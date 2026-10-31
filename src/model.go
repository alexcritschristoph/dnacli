package main

import (
	"fmt"
	"os"
	"os/exec"
	"sort"
	"strconv"
	"strings"
	"time"

	osc52 "github.com/aymanbagabas/go-osc52/v2"
	"github.com/biogo/hts/bam"
	tea "github.com/charmbracelet/bubbletea"
)

type viewState struct {
	start int // leftmost reference position, zero-based
	bpp   int // bases per column
}

type model struct {
	src            *Source
	tabs           []*Source // every alignment file opened, in the order they were opened
	tab            int       // which of them the view is showing
	ref            *Reference
	refPath        string
	variants       []Variant
	filter         Filter
	highlight      int  // the SNP highlighting chosen, -1 for none
	color          int  // the color scheme chosen, -1 for none
	palette        int  // the set of colors the schemes draw from
	schemeRow      int  // the coloring panel row the arrow keys are on
	coloringLocked bool // whether the panel has stopped previewing as the cursor moves

	width, height int
	frame         grid // the screen we draw into, kept between frames
	exporting     bool // drawing a picture, so overlays stay out of it
	contig        int
	view          viewState
	shownView     viewState // the view the screen on display describes
	shownSpan     int       // bases the screen on display covers
	rowOffset     int
	wheelTurned   time.Time // when the wheel last turned, which tells a mouse from a track pad
	wheelActed    time.Time // when a turn of the wheel last did something
	movePan       int       // bases the keys have asked for, not yet shown
	moveZoom      int       // ladder steps the keys have asked for, not yet shown
	moveAnchor    int       // the column a pending zoom holds still
	moveScheduled bool      // a timer is booked to show the pending movement
	async         bool      // whether a movement reads the file off the drawing thread
	loading       bool      // a screen is being built while the last one stays up
	stale         bool      // the model changed while a screen was being built
	quiet         bool      // pan and zoom change the view without reading the file
	spin          int       // which frame of the loading mark is shown

	reads      []*Alignment
	held       window       // the reads kept around the view, so a pan only asks for the new strip
	cols       []columnSpec // what each screen column shows
	inserts    map[int]int  // columns given over to insertions, by the position they follow
	readLength int          // typical read length here, which decides how far we can zoom out
	rowCount   int
	readsSeen  int // reads the region held, before the filters
	truncated  bool
	sites      map[int]site
	consensus  []byte // the base most reads carry, by position from the left edge
	counted    counts // the per-base tallies of the last stretch counted
	cov        []int
	covPeak    int

	panel        string // "", variants, contigs, colors, filters, help
	panelFocus   bool   // whether the arrow keys belong to the panel or the reads
	variantIndex int
	contigIndex  int
	filterIndex  int
	command      string // what has been typed into the bar along the top
	choice       int    // the highlighted line of its suggestions

	hoverPos     int
	marked       int // where the go command left the mark, or -1
	hover        *Alignment
	impacts      map[int]Impact
	mouseOn      bool
	allFeatures  bool // whether the gene row shows every GenBank feature, not just genes and CDSs
	message      string
	messageUntil time.Time
	flash        *Alignment // a read lit up for a moment, to show what was copied
	flashUntil   time.Time
	bars         []bar    // the labelled strips below the pileup
	tabHits      []tabHit // where the tabs were last drawn, so a click can find them
	pending      tea.Cmd  // work a command leaves for the screen to start
}

// tabHit is the span of one drawn tab along the top bar.
type tabHit struct {
	from, to int
	tab      int
}

// profileMsg carries a file's finished coverage sketches back to the screen.
type profileMsg struct {
	source   *Source
	profiles map[string][]float64
}

// profileBins is the fewest coverage bins a contig gets. Longer contigs get
// about one bin every 200 bases, so the depth track stays smooth out to the
// zoom where the reads themselves stop being drawn.
const profileBins = 2048

func profileBinsFor(length int) int {
	return clamp(length/200, profileBins, 1<<16)
}

var zoomLadder = []int{1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000}

func newModel(src *Source, ref *Reference, refPath string, variants []Variant) *model {
	m := &model{src: src, tabs: []*Source{src}, ref: ref, refPath: refPath, variants: variants,
		filter: defaultFilter(), view: viewState{0, 1}, mouseOn: true, variantIndex: -1,
		impacts: map[int]Impact{}, readLength: 150, width: 100, height: 30, hoverPos: -1, marked: -1}
	chosenPalette = palettes[0]
	for i, c := range src.Contigs { // open on the first contig that has reads
		if c.Reads > 0 {
			m.contig, m.contigIndex = i, i
			break
		}
	}
	m.highlight = namedIndex(highlights, "Mismatch (from consensus)") // without a reference, the reads are the guide
	if ref != nil {
		m.highlight = namedIndex(highlights, "Mismatch (from reference)")
	}
	m.schemeRow, m.color = m.highlight, -1
	if len(variants) > 0 {
		m.panel = "variants"
		m.measureImpacts()
	} else if len(src.Contigs) > 1 {
		m.panel = "contigs"
	}
	if ref != nil && ref.UnsupportedTable != 0 {
		m.notify("Genetic code table %d is not built in; the standard code is used instead.", ref.UnsupportedTable)
	}
	return m
}

func (m *model) contigs() []Contig { return m.src.Contigs }

// measureImpacts works out what every variant does to the reference's proteins.
// The CDS positions behind each call are cached, so asking about thousands of
// variants in the same genes does not rebuild the same list every time.
func (m *model) measureImpacts() {
	for i := range m.variants {
		v := &m.variants[i]
		v.Impact = m.ref.CodingImpact(v.Contig, v.Start, v.Ref, v.Alt)
	}
}

// Init measures coverage across every contig that has reads, in the background.
func (m *model) Init() tea.Cmd {
	var work []tea.Cmd
	for _, src := range m.tabs {
		work = append(work, sketchCoverage(src))
	}
	work = append(work, m.fade()) // a message set before the screen opened still fades
	return tea.Sequence(work...)
}

// sketchCoverage measures the coverage of every contig in a file. It opens the
// file once, in the background, so that the screen keeps its own handle and a
// file of a thousand contigs is not opened a thousand times.
func sketchCoverage(src *Source) tea.Cmd {
	return func() tea.Msg {
		own, err := OpenAlignments(src.Path)
		if err != nil {
			return nil
		}
		defer own.Close()
		if own.reader != nil {
			own.reader.Omit(bam.AllVariableLengthData) // a sketch needs positions, not bases
		}
		profiles := map[string][]float64{}
		for _, c := range own.Contigs {
			if c.Reads > 0 {
				profiles[c.Name] = own.CoverageProfile(c.Name, profileBinsFor(c.Length))
			}
		}
		return profileMsg{src, profiles}
	}
}

func (m *model) Update(msg tea.Msg) (tea.Model, tea.Cmd) {
	switch msg := msg.(type) {
	case tea.WindowSizeMsg:
		m.width, m.height = msg.Width, msg.Height
		m.reload()
	case profileMsg:
		for _, s := range m.tabs {
			if s != msg.source {
				continue
			}
			for i, c := range s.Contigs {
				if p := msg.profiles[c.Name]; p != nil {
					s.Contigs[i].Profile = p
				}
			}
		}
		if msg.source == m.src {
			m.reload()
		}
	case tea.MouseMsg:
		if m.loading {
			return m, nil // the frame on screen is frozen, so there is nothing new to point at
		}
		return m, m.handleMouse(msg)
	case tea.KeyMsg:
		return m, tea.Batch(m.handleKey(msg), m.fade())
	case fadeMsg:
		if !time.Now().Before(m.flashUntil) {
			m.flash = nil
		}
		if !time.Now().Before(m.messageUntil) {
			m.message = ""
		}
		return m, m.fade()
	case moveMsg:
		return m, m.showMove()
	case loadedMsg:
		if m.stale { // something other than a movement changed while the screen was built
			m.held = msg.model.held // keep the reads the file gave back
			m.loading, m.stale = false, false
			m.reload()
			return m, nil
		}
		// Keep what the hand did while the screen was built.
		pan, zoom, anchor, scheduled := m.movePan, m.moveZoom, m.moveAnchor, m.moveScheduled
		offset, message, messageUntil := m.rowOffset, m.message, m.messageUntil
		flash, flashUntil := m.flash, m.flashUntil
		*m = *msg.model
		m.movePan, m.moveZoom, m.moveAnchor, m.moveScheduled = pan, zoom, anchor, scheduled
		m.rowOffset = clamp(offset, 0, max(0, m.rowCount-m.layout().pileupRows))
		m.message, m.messageUntil, m.flash, m.flashUntil = message, messageUntil, flash, flashUntil
		if m.movePan != 0 || m.moveZoom != 0 { // keys arrived while the screen was built
			return m, tea.Tick(movePause, func(time.Time) tea.Msg { return moveMsg{} })
		}
		return m, nil
	case spinMsg:
		if m.loading {
			m.spin++
			return m, m.spinTick()
		}
		return m, nil
	}
	return m, nil
}

func (m *model) handleKey(msg tea.KeyMsg) tea.Cmd {
	key := msg.String()
	if key == "ctrl+c" || (key == "q" && m.command == "") {
		return tea.Quit
	}
	if m.command != "" {
		if m.editCommand(msg) {
			return m.takePending()
		}
		m.command, m.choice = "", 0 // any other key leaves the bar and works the view
	}
	step := max(1, m.span()*7/50) // a little under a fifth of the screen
	switch key {
	case "t":
		m.nextTab()
	case "esc":
		m.panel, m.panelFocus = "", false
		m.reload()
	case "?":
		m.togglePanel("help")
	case "left":
		return m.moveSideways(-1, -step)
	case "right":
		return m.moveSideways(1, step)
	case "up":
		m.moveUpDown(-1)
	case "down":
		m.moveUpDown(1)
	case "pgup":
		m.scroll(-m.layout().pileupRows)
	case "pgdown":
		m.scroll(m.layout().pileupRows)
	case "home":
		m.goTo(0)
	case "end":
		m.goTo(m.contigs()[m.contig].Length)
	case "+", "=":
		return m.zoomLater(-1, m.layout().columns/2)
	case "-", "_":
		return m.zoomLater(1, m.layout().columns/2)
	case "0":
		m.view.bpp = m.fitZoom()
		m.goTo(m.contigs()[m.contig].Length / 2)
	case "z":
		m.view.bpp = 1
		m.reload()
	case "a":
		m.allFeatures = !m.allFeatures
	case "m":
		m.mouseOn = !m.mouseOn
		if m.mouseOn {
			m.notify("Mouse on: hover on reads and click them to copy")
			return tea.EnableMouseAllMotion
		}
		m.notify("Mouse off: select and copy text.")
		m.clearHover()
		return tea.DisableMouse
	case "enter":
		if m.panelFocus && m.panel == "colors" {
			m.pressColoring()
			return nil
		}
		m.panelFocus = m.panel != "" && m.panel != "help"
	default:
		if msg.Type == tea.KeyRunes || msg.Type == tea.KeySpace { // anything else starts a command
			m.editCommand(msg)
		}
	}
	return nil
}

// moveUpDown is the up and down keys: they work the open panel when it has the
// focus, and scroll the reads when it does not.
func (m *model) moveUpDown(delta int) {
	if !m.panelFocus || m.panel == "" {
		m.scroll(delta)
		return
	}
	switch m.panel {
	case "variants":
		m.showVariant(clamp(m.variantIndex+delta, 0, len(m.variants)-1))
	case "contigs":
		m.selectContig(m.contigIndex + delta)
	case "colors":
		m.focusColoring(m.schemeRow + delta)
	case "filters":
		m.filterIndex = clamp(m.filterIndex+delta, 0, len(filterRows)-1)
	}
}

// moveSideways is the left and right keys: they change the chosen filter when
// the filter panel has the focus, and pan otherwise.
func (m *model) moveSideways(delta, bases int) tea.Cmd {
	if m.panelFocus && m.panel == "filters" {
		filterRows[m.filterIndex].move(&m.filter, delta)
		m.reload()
		return nil
	}
	return m.panLater(bases)
}

// openPanel shows a panel and hands it the arrow keys.
func (m *model) openPanel(name string) {
	m.panel, m.panelFocus = name, true
	if name == "colors" {
		m.coloringLocked = false // a freshly opened panel previews as the cursor moves
	}
	m.reload()
	if m.layout().panel == 0 {
		m.notify("The window is too narrow for the side panel.")
	}
}

func (m *model) togglePanel(name string) {
	if m.panel == name {
		m.panel = ""
	} else {
		m.panel = name
		if name == "colors" {
			m.coloringLocked = false
		}
	}
	m.reload()
}

// notify writes a line above the details, which fades on its own.
func (m *model) notify(format string, args ...any) {
	m.message, m.messageUntil = fmt.Sprintf(format, args...), time.Now().Add(messageLife)
}

// fadeMsg arrives once a message or a flash has had its time on screen.
type fadeMsg struct{}

const (
	messageLife = 5 * time.Second
	flashLife   = 400 * time.Millisecond
)

// fade asks to be woken when whatever is on screen next runs out of time.
func (m *model) fade() tea.Cmd {
	var until time.Time
	if m.flash != nil {
		until = m.flashUntil
	}
	if m.message != "" && (until.IsZero() || m.messageUntil.Before(until)) {
		until = m.messageUntil
	}
	if until.IsZero() {
		return nil
	}
	return tea.Tick(time.Until(until)+time.Millisecond, func(time.Time) tea.Msg { return fadeMsg{} })
}

// moving around ------------------------------------------------------------

func (m *model) span() int { return m.layout().columns * m.view.bpp }

func (m *model) pan(bases int) {
	m.view.start += bases
	m.reload()
}

// A held arrow key sends a key every 20 to 30 ms, faster than a screen of reads
// can be loaded. So a pan or a zoom only counts itself up, and a timer a moment
// later shows what is pending in one go: the keys drain at once instead of the
// view replaying them after the hand has stopped, and a pan that is followed by a
// zoom or a pan the other way never happens at all.
const movePause = 20 * time.Millisecond

type moveMsg struct{}

// panLater counts up a pan for the timer to show.
func (m *model) panLater(bases int) tea.Cmd {
	m.movePan += bases
	m.moveZoom = 0 // the newest movement is the one that speaks
	return m.scheduleMove()
}

// zoomLater counts up a zoom for the timer to show, holding still the base under
// the anchor column.
func (m *model) zoomLater(steps, anchor int) tea.Cmd {
	m.moveZoom += steps
	m.movePan = 0
	m.moveAnchor = anchor
	return m.scheduleMove()
}

// scheduleMove books the one timer that shows the pending movement.
func (m *model) scheduleMove() tea.Cmd {
	if m.moveScheduled {
		return nil
	}
	m.moveScheduled = true
	return tea.Tick(movePause, func(time.Time) tea.Msg { return moveMsg{} })
}

// showMove pans and zooms as far as the keys have asked since the last frame.
// The view moves only when a screen for it is ready, so the keys never leave the
// picture describing a stretch it no longer shows.
func (m *model) showMove() tea.Cmd {
	m.moveScheduled = false
	if m.movePan == 0 && m.moveZoom == 0 {
		return nil // nothing moved, so nothing to build
	}
	if m.loading {
		return nil // the keys keep counting up; the screen on its way will be followed by this
	}
	pan, zoom, anchor := m.movePan, m.moveZoom, m.moveAnchor
	m.movePan, m.moveZoom = 0, 0
	m.quiet = true // the view moves now; the screen is built once below
	if pan != 0 {
		m.pan(pan)
	}
	if zoom != 0 {
		m.zoom(zoom, anchor)
	}
	m.quiet = false
	if !m.async {
		m.reload()
		return nil
	}
	target := m.view
	m.view = m.shownView // the frame on screen stays the one the model describes
	m.loading, m.stale = true, false
	clone := *m
	clone.view = target
	return tea.Batch(func() tea.Msg {
		clone.loading = false
		clone.reload()
		return loadedMsg{view: target, model: &clone}
	}, m.spinTick())
}

// loadedMsg carries a screen built off the drawing thread.
type loadedMsg struct {
	view  viewState
	model *model
}

// spinMsg asks for the next frame of the loading mark.
type spinMsg struct{}

func (m *model) spinTick() tea.Cmd {
	return tea.Tick(90*time.Millisecond, func(time.Time) tea.Msg { return spinMsg{} })
}

// showRange frames a stretch of the reference, with a little room either side.
func (m *model) showRange(start, end int) {
	columns := m.layout().columns
	m.view.bpp = clamp(((end-start)*11/10+columns-1)/columns, 1, m.fitZoom())
	m.goTo((start + end) / 2)
}

// goTo centres the view on a position.
func (m *model) goTo(pos int) {
	m.view.start = pos - m.span()/2
	m.reload()
}

func (m *model) scroll(rows int) {
	m.rowOffset = clamp(m.rowOffset+rows, 0, max(0, m.rowCount-m.layout().pileupRows))
}

// A mouse sends one turn of the wheel for each click of it, a track pad tens of
// turns a second, so acting on every one is either too slow or far too fast.
// Turns further apart than mouseClick come from a hand and all act; a track pad's
// stream of them acts at most once every wheelPause.
const mouseClick, wheelPause = 25 * time.Millisecond, 80 * time.Millisecond

// zoom moves one step along the ladder, holding still whatever the anchor
// column is showing, so that zooming with the wheel closes in on the pointer.
// A zoom clamped to the whole contig sits between rungs, so it steps to the
// neighbouring rung rather than skipping the one it landed past.
func (m *model) zoom(steps, anchor int) {
	bpp := m.view.bpp
	for ; steps > 0; steps-- {
		bpp = nextZoom(bpp, 1, m.fitZoom())
	}
	for ; steps < 0; steps++ {
		bpp = nextZoom(bpp, -1, m.fitZoom())
	}
	if bpp == m.view.bpp {
		return // already at that end of the ladder; a trackpad's tail must not pan
	}
	pos := m.view.start + anchor*m.view.bpp
	if anchor < len(m.cols) {
		pos = m.cols[anchor].pos
	}
	m.view.bpp = bpp
	m.view.start = pos - anchor*bpp
	m.reload()
}

// nextZoom is the rung of the ladder one step away from bpp, never past the
// whole-contig zoom.
func nextZoom(bpp, steps, fit int) int {
	i := sort.SearchInts(zoomLadder, bpp)
	if steps < 0 {
		if i == 0 {
			return bpp
		}
		return zoomLadder[i-1]
	}
	if i < len(zoomLadder) && zoomLadder[i] == bpp {
		i++
	}
	if i == len(zoomLadder) {
		return bpp
	}
	return min(zoomLadder[i], fit)
}

// fitZoom is the widest useful zoom: the whole contig on one screen.
func (m *model) fitZoom() int {
	columns := m.layout().columns
	return max(1, (m.contigs()[m.contig].Length+columns-1)/columns)
}

// takePending hands the screen work a command left behind, once.
func (m *model) takePending() tea.Cmd {
	work := m.pending
	m.pending = nil
	return work
}

// nextTab moves on to the next alignment file, around to the first after the
// last.
func (m *model) nextTab() {
	if len(m.tabs) > 1 {
		m.switchTab((m.tab + 1) % len(m.tabs))
	}
}

// switchTab shows another alignment file. The view keeps its position when the
// new file has the same contig, and otherwise opens on the first contig with
// reads there.
func (m *model) switchTab(tab int) {
	if tab = clamp(tab, 0, len(m.tabs)-1); tab == m.tab {
		return
	}
	name := m.contigs()[m.contig].Name
	m.tab, m.src = tab, m.tabs[tab]
	m.contig, m.contigIndex = 0, 0
	for i, c := range m.contigs() {
		if c.Name == name {
			m.contig, m.contigIndex = i, i
			break
		}
	}
	if m.contigs()[m.contig].Name != name {
		for i, c := range m.contigs() {
			if c.Reads > 0 {
				m.contig, m.contigIndex = i, i
				break
			}
		}
	}
	m.clearHover()
	m.marked = -1
	m.reload()
}

func (m *model) selectContig(i int) {
	if i = clamp(i, 0, len(m.contigs())-1); i == m.contig {
		return
	}
	m.contig, m.contigIndex = i, i
	m.view, m.marked = viewState{0, 1}, -1
	m.rowOffset = 0
	m.reload()
}

func (m *model) showVariant(i int) {
	if i < 0 || i >= len(m.variants) {
		return
	}
	m.variantIndex, m.marked = i, -1 // the variant takes the mark
	v := m.variants[i]
	for c, contig := range m.contigs() {
		if contig.Name == v.Contig && c != m.contig {
			m.selectContig(c)
		}
	}
	m.goTo(v.Start)
}

func (m *model) selectedVariant() *Variant {
	if m.variantIndex < 0 || m.variantIndex >= len(m.variants) {
		return nil
	}
	return &m.variants[m.variantIndex]
}

// parseGoto accepts "1200", "contig", "contig:1200" and "contig:1200-1400".
func (m *model) parseGoto(text string) {
	text = strings.TrimSpace(text)
	name, where, hasContig := strings.Cut(text, ":")
	if hasContig {
		found := false
		for i, c := range m.contigs() {
			if c.Name == name {
				m.selectContig(i)
				found = true
			}
		}
		if !found {
			m.notify("No contig called %q.", name)
			return
		}
	} else {
		where = text
	}
	from, to, isRange := strings.Cut(where, "-")
	start, err := strconv.Atoi(strings.ReplaceAll(strings.TrimSpace(from), ",", ""))
	if err != nil {
		if !hasContig {
			m.notify("Type a position, or contig:position.")
		}
		return
	}
	if end, err := strconv.Atoi(strings.ReplaceAll(strings.TrimSpace(to), ",", "")); isRange && err == nil && end > start {
		m.view.bpp = clamp((end-start+m.layout().columns-1)/m.layout().columns, 1, m.fitZoom())
		m.goTo((start + end) / 2)
		return
	}
	m.marked = start - 1 // going to a spot marks it, the way choosing a variant does
	m.goTo(m.marked)
}

// clipboardTools reach the clipboard of the machine dnacli is running on,
// tried in order. Over ssh that is the wrong machine, so they are skipped.
var clipboardTools = [][]string{
	{"wl-copy"},
	{"xclip", "-selection", "clipboard"},
	{"xsel", "--input", "--clipboard"},
	{"pbcopy"},
	{"clip.exe"},
}

// copyText puts text on the clipboard by every route open to it, and says
// whether one of them is known to have worked. The escape sequence goes out
// whatever happens: it is the only route that reaches the terminal you are
// sitting at when dnacli is running somewhere else.
func copyText(text string) bool {
	remote := os.Getenv("SSH_TTY") != ""
	trusted := false
	if os.Getenv("TMUX") != "" {
		trusted = pipeTo(text, "tmux", "load-buffer", "-w", "-")
	}
	if !trusted && !remote {
		for _, tool := range clipboardTools {
			if trusted = pipeTo(text, tool[0], tool[1:]...); trusted {
				break
			}
		}
	}
	osc52.New(text).WriteTo(os.Stderr)
	return trusted || remote // over ssh the escape is the only way out, and usually works
}

func pipeTo(text, name string, args ...string) bool {
	path, err := exec.LookPath(name)
	if err != nil {
		return false
	}
	cmd := exec.Command(path, args...)
	cmd.Stdin = strings.NewReader(text)
	return cmd.Run() == nil
}

// mouse --------------------------------------------------------------------

func (m *model) clearHover() { m.hover, m.hoverPos = nil, -1 }

func (m *model) handleMouse(msg tea.MouseMsg) tea.Cmd {
	l := m.layout()
	if msg.Action == tea.MouseActionPress {
		m.command, m.choice = "", 0 // a click anywhere, or a turn of the wheel, leaves the bar
	}
	if msg.Action == tea.MouseActionPress && msg.Button == tea.MouseButtonLeft && msg.Y == 0 {
		for _, hit := range m.tabHits {
			if msg.X >= hit.from && msg.X <= hit.to {
				m.switchTab(hit.tab)
				return nil
			}
		}
	}
	if msg.Action == tea.MouseActionPress && msg.Button == tea.MouseButtonLeft &&
		msg.Y == trackTop+l.geneRows && msg.X == gutterWidth-3 && l.geneRows > 0 &&
		(m.allFeatures || m.hiddenFeatures(l.columns)) {
		m.allFeatures = !m.allFeatures
		return nil
	}
	if l.panel > 0 && msg.X >= m.width-l.panel {
		m.clearHover() // the pointer has left the reads, so the line falls back to the variant
		if msg.Action == tea.MouseActionPress && msg.Button == tea.MouseButtonLeft {
			m.panelFocus = true // clicking the panel gives it the arrow keys
			m.clickPanel(msg.Y, l)
		}
		return nil
	}
	switch msg.Button {
	case tea.MouseButtonWheelUp, tea.MouseButtonWheelDown:
		stream := time.Since(m.wheelTurned) < mouseClick
		m.wheelTurned = time.Now()
		if stream && time.Since(m.wheelActed) < wheelPause {
			return nil
		}
		m.wheelActed = time.Now()
		turn := map[bool]int{true: -1, false: 1}[msg.Button == tea.MouseButtonWheelUp]
		switch {
		case msg.Ctrl:
			m.scroll(turn * 3)
		case msg.Shift:
			return m.panLater(turn * max(1, m.span()/8))
		default:
			return m.zoomLater(turn, clamp(msg.X-gutterWidth, 0, l.columns-1))
		}
		return nil
	}
	defer m.refreshBars() // the strips below the pileup follow the pointer
	read := m.readAt(msg.X, msg.Y, l)
	if msg.Action == tea.MouseActionMotion || msg.Action == tea.MouseActionRelease {
		m.hover, m.hoverPos = read, m.positionAt(msg.X)
		return nil
	}
	if msg.Action != tea.MouseActionPress || msg.Button != tea.MouseButtonLeft {
		return nil
	}
	m.panelFocus = false // clicking the reads takes the arrow keys back
	if read == nil {
		return nil
	}
	return m.copyRead(read)
}

// copyRead puts a read on the clipboard as FASTA, and says so on screen: the
// clipboard itself gives no sign that anything happened.
func (m *model) copyRead(a *Alignment) tea.Cmd {
	if copyText(">" + a.Name + "\n" + a.Sequence) {
		m.notify("read %s copied to clipboard", a.Name)
	} else {
		m.notify("read %s copied — this terminal may not accept it", a.Name)
	}
	m.flash, m.flashUntil = a, time.Now().Add(flashLife)
	return m.fade()
}

// positionAt is the reference position a screen column shows.
func (m *model) positionAt(x int) int {
	if col := x - gutterWidth; col >= 0 && col < len(m.cols) {
		return m.cols[col].pos
	}
	return -1
}

func (m *model) readAt(x, y int, l layout) *Alignment {
	row := y - l.pileupTop + m.rowOffset
	col := x - gutterWidth
	if y < l.pileupTop || y >= l.pileupTop+l.pileupRows || col < 0 || col >= l.columns {
		return nil
	}
	for _, a := range m.reads {
		if a.Row == row && col >= m.column(a.Start) && col <= m.column(a.End-1) {
			return a
		}
	}
	return nil
}

func (m *model) clickPanel(y int, l layout) {
	row := y - l.listTop
	switch m.panel {
	case "variants":
		offset := max(0, min(m.variantIndex-l.listRows/2, len(m.variants)-l.listRows))
		m.showVariant(offset + row)
	case "contigs":
		visible := l.listRows / 3
		offset := max(0, min(m.contigIndex-visible/2, len(m.contigs())-visible))
		m.selectContig(offset + row/3)
	case "colors":
		if i := m.coloringOffset(l) + row; i >= 0 && i < len(coloringPanel) && coloringPanel[i].entry >= 0 {
			m.pickColoring(coloringPanel[i].entry)
		}
	case "filters":
		if i := row / 2; i >= 0 && i < len(filterRows) {
			m.filterIndex = i
		}
	}
}

// window is the stretch of reference whose reads are held, so that a pan asks
// the file only for the strip that has scrolled into view.
type window struct {
	src        *Source
	start, end int
	contig     int
	detail     bool
	filter     Filter
	reads      []*Alignment
	seen       int
	truncated  bool
}

// windowReads returns the reads overlapping a stretch, reading only the strips
// the window does not already hold. It keeps a screen either side of the view,
// so a small pan costs nothing and a large jump starts over.
func (m *model) windowReads(contig, start, end int, detail bool) ([]*Alignment, int, bool) {
	w := &m.held
	name := m.contigs()[contig].Name
	// A window with a base for every position can answer a wider view that only
	// wants the differences, so zooming out does not have to read it all again.
	stale := detail && !w.detail
	if w.src != m.src || w.contig != contig || stale || w.filter != m.filter ||
		w.start >= w.end || end <= w.start || start >= w.end {
		*w = window{src: m.src, contig: contig, detail: detail, filter: m.filter, start: start, end: end}
		w.reads, w.seen, w.truncated = m.src.Reads(name, start, end, m.filter, m.ref, detail)
		m.counted = counts{} // the old tallies belong to another stretch
		return w.reads, w.seen, w.truncated
	}
	left, right, added := w.start, w.end, false
	share := func(width int) int { return max(256, maxAlignments*width/max(1, end-start)) }
	// When a pan has only exposed a sliver, read a whole screen past the edge so
	// the next few pans land inside what is held. A zoom exposes much more at
	// once, so there it reads just what the view needs.
	tile := m.span()
	if end-start != m.shownSpan {
		tile = 0 // a zoom or a resize takes much more than a pan sliver
	}
	if start < left { // panned left: read what ended before the held stretch
		from := start
		if left-start < tile/2 {
			from = max(0, left-tile)
		}
		reads, seen, truncated := m.src.readsUpTo(name, from, left, share(left-from), m.filter, m.ref, detail)
		for _, a := range reads {
			if a.End <= left {
				w.reads = append(w.reads, a)
			}
		}
		w.start, w.seen, w.truncated, added = from, w.seen+seen, w.truncated || truncated, true
	}
	if end > right { // panned right: read what starts after the held stretch
		to := end
		if end-right < tile/2 {
			to = min(m.contigs()[contig].Length, right+tile)
		}
		reads, seen, truncated := m.src.readsUpTo(name, right, to, share(to-right), m.filter, m.ref, detail)
		for _, a := range reads {
			if a.Start >= right {
				w.reads = append(w.reads, a)
			}
		}
		w.end, w.seen, w.truncated, added = to, w.seen+seen, w.truncated || truncated, true
	}
	if added {
		sort.SliceStable(w.reads, func(i, j int) bool { return w.reads[i].Start < w.reads[j].Start })
		w.detail = w.detail && detail // the strips just read may lack bases
	}
	// Let go of anything more than a screen away from the view.
	lo, hi := max(0, start-m.span()), min(m.contigs()[contig].Length, end+m.span())
	if lo > w.start || hi < w.end {
		kept := w.reads[:0]
		for _, a := range w.reads {
			if a.End > lo && a.Start < hi {
				kept = append(kept, a)
			}
		}
		w.reads, w.start, w.end = kept, max(w.start, lo), min(w.end, hi)
	}
	out := make([]*Alignment, 0, len(w.reads))
	for _, a := range w.reads {
		if a.End > start && a.Start < end {
			out = append(out, a)
		}
	}
	return out, w.seen, w.truncated
}

// counts holds the per-base tallies of the stretch they were worked out for, so
// that a pan only counts the columns that have scrolled into view.
type counts struct {
	src        *Source
	contig     int
	filter     Filter
	start, end int
	counts     [][4]int
	depth      []int
}

// countBases tallies the bases of a stretch, reusing the columns it already has
// tallies for and counting only the rest.
func (m *model) countBases(reads []*Alignment, start, end int) ([][4]int, []int) {
	c := &m.counted
	if c.src != m.src || c.contig != m.contig || c.filter != m.filter {
		*c = counts{src: m.src, contig: m.contig, filter: m.filter}
	}
	tallies, depth := make([][4]int, end-start), make([]int, end-start)
	if from, to := max(start, c.start), min(end, c.end); from < to {
		copy(tallies[from-start:], c.counts[from-c.start:to-c.start])
		copy(depth[from-start:], c.depth[from-c.start:to-c.start])
	}
	if start < c.start {
		countRange(tallies, depth, start, reads, start, min(end, c.start))
	}
	if end > c.end {
		countRange(tallies, depth, start, reads, max(start, c.end), end)
	}
	c.start, c.end, c.counts, c.depth = start, end, tallies, depth
	return tallies, depth
}

// loading ------------------------------------------------------------------

func (m *model) reload() {
	if m.quiet {
		return // the movement below will build one screen for all of it
	}
	if m.loading {
		m.stale = true // the screen being built is already out of date
		return
	}
	m.movePan, m.moveZoom, m.moveScheduled = 0, 0, false // whatever was pending is now out of date
	l := m.layout()
	contig := m.contigs()[m.contig]
	m.view.bpp = clamp(m.view.bpp, 1, m.fitZoom())
	span := l.columns * m.view.bpp
	m.view.start = clamp(m.view.start, 0, max(0, contig.Length-span)) // hard stop at both ends
	end := min(m.view.start+span, contig.Length)

	if !m.showReads() {
		// Too far out for single reads: draw the depth sketch and nothing else.
		m.reads, m.rowCount, m.readsSeen, m.truncated = nil, 0, 0, false
		m.buildColumns(l.columns)
		m.sites = map[int]site{}
		m.consensus = nil
		m.cov, m.covPeak = depthFromProfile(contig.Profile, contig.Length, m.view.start, m.view.bpp, l.columns)
		m.impacts = map[int]Impact{}
		m.shownView, m.shownSpan = m.view, span
		return
	}
	// A screen wider than baseDetailSpan cannot show single bases, so keep only
	// what it draws: read spans, and the bases that differ from the reference.
	// The consensus of the reads is a per-base idea too, so it waits for a
	// closer zoom rather than costing a pass over every base here.
	detail := span <= baseDetailSpan
	m.reads, m.readsSeen, m.truncated = m.windowReads(m.contig, m.view.start, end, detail)
	m.measureReads()
	m.buildColumns(l.columns)
	m.rowCount = packRows(m.reads, max(1, m.view.bpp), span)
	if detail {
		counts, depth := m.countBases(m.reads, m.view.start, end)
		m.consensus = consensusBases(counts)
		m.sites = callSites(counts, depth, m.ref, contig.Name, m.view.start)
	} else {
		m.consensus, m.sites = nil, map[int]site{} // nothing this wide is drawn from single bases
		m.counted = counts{}
	}
	m.cov, m.covPeak = coverage(m.reads, m.cols, m.view.bpp)
	m.rowOffset = clamp(m.rowOffset, 0, max(0, m.rowCount-l.pileupRows))
	m.clearHover()
	m.impacts = map[int]Impact{}
	m.shownView, m.shownSpan = m.view, span
}

// baseColor is the color of one drawn base: the scheme's own color when the
// read shows what belongs at that position, and something louder when it does
// not. Expected means the reference base, or the base most reads carry when
// the coloring asks for the consensus.
func (m *model) baseColor(a *Alignment, pos int, base, refBase byte) (string, bool) {
	if base == '-' || base == ' ' {
		return readColor(m.color, a, pos, m.readLength), false
	}
	expected := m.expectedBase(pos, refBase)
	if expected == 0 || base == expected {
		return readColor(m.color, a, pos, m.readLength), false
	}
	if m.color >= 0 && colors[m.color] == "Protein change" && m.ref != nil {
		return impactColor(m.impactAt(pos, refBase, base)), true
	}
	return chosenPalette.bases[base], true
}

// expectedBase is the base a position should hold for the SNP highlighting in
// force: the reference, or the consensus of the reads. Zero means there is
// nothing to compare against, or no highlighting is on.
func (m *model) expectedBase(pos int, refBase byte) byte {
	switch {
	case m.highlight < 0:
		return 0
	case highlights[m.highlight] == "Mismatch (from consensus)":
		return m.consensusBaseAt(pos)
	}
	return refBase
}

// fromConsensus says whether the SNP highlighting compares against the
// consensus of the reads rather than the reference.
func (m *model) fromConsensus() bool {
	return m.highlight >= 0 && highlights[m.highlight] == "Mismatch (from consensus)"
}

// consensusBaseAt is the base most reads carry at a position, or 0 if there is
// no consensus there.
func (m *model) consensusBaseAt(pos int) byte {
	if i := pos - m.view.start; i >= 0 && i < len(m.consensus) {
		return m.consensus[i]
	}
	return 0
}

// focusColoring moves the arrow keys to another entry of the coloring panel.
// Until a choice has been locked, the entry under the cursor is also the one in
// force, so that moving through the panel previews each setting. Once locked,
// only the cursor moves.
func (m *model) focusColoring(entry int) {
	m.schemeRow = clamp(entry, 0, len(highlights)+len(colors)+len(palettes)-1)
	if m.coloringLocked {
		return
	}
	switch {
	case m.schemeRow < len(highlights):
		m.setHighlight(m.schemeRow)
	case m.schemeRow < len(highlights)+len(colors):
		m.color = m.schemeRow - len(highlights)
	default:
		m.choosePalette(m.schemeRow - len(highlights) - len(colors))
	}
}

// pressColoring is what enter does in the coloring panel. The first press locks
// the highlighted entry, after which the arrows only move the cursor. A press
// on the locked entry clears it; a press on another entry locks that one
// instead.
func (m *model) pressColoring() {
	switch {
	case m.schemeRow < len(highlights):
		if m.coloringLocked && m.highlight == m.schemeRow {
			m.setHighlight(-1)
			return
		}
		m.setHighlight(m.schemeRow)
	case m.schemeRow < len(highlights)+len(colors):
		i := m.schemeRow - len(highlights)
		if m.coloringLocked && m.color == i {
			m.color = -1
			return
		}
		m.color = i
	default: // a palette is always in force, so pressing it again keeps it
		m.choosePalette(m.schemeRow - len(highlights) - len(colors))
	}
	m.coloringLocked = true
}

// choosePalette repaints every scheme with another set of colors, leaving the
// chosen highlighting and color scheme as they are.
func (m *model) choosePalette(i int) {
	m.palette = i
	chosenPalette = palettes[i]
}

// pickColoring moves the cursor to the entry that was clicked, and presses it.
func (m *model) pickColoring(entry int) {
	m.schemeRow = entry
	m.pressColoring()
}

// coloringOffset is the first line of the coloring panel on screen: the least
// scrolling that keeps the cursor in view.
func (m *model) coloringOffset(l layout) int {
	line := 0
	for i, row := range coloringPanel {
		if row.entry == m.schemeRow {
			line = i
		}
	}
	return clamp(line-l.listRows+1, 0, max(0, len(coloringPanel)-l.listRows))
}

// setHighlight chooses the SNP highlighting, or none, and works out the marks
// it needs.
func (m *model) setHighlight(i int) {
	m.highlight = i
}

// chosen says whether an entry of the coloring panel is the one in force.
func (m *model) chosen(entry int) bool {
	if entry < len(highlights) {
		return entry == m.highlight
	}
	entry -= len(highlights)
	if entry < len(colors) {
		return entry == m.color
	}
	return entry-len(colors) == m.palette
}

// coloringSummary names what the reads are drawn with: whichever parts of the
// coloring panel are on, and the palette when it is not the default.
func (m *model) coloringSummary() string {
	var summary string
	switch {
	case m.color < 0 && m.highlight < 0:
		summary = "Color by: none"
	case m.color < 0:
		summary = "Color by: " + highlights[m.highlight]
	case m.highlight < 0:
		summary = "Color by: " + colors[m.color]
	default:
		summary = "Color by: " + colors[m.color] + " + " + highlights[m.highlight]
	}
	if m.palette != 0 {
		summary += " · " + palettes[m.palette].name + " palette"
	}
	return summary
}

// impactAt works out what one base change does to the protein, remembering
// answers because a screen full of reads asks the same question many times.
func (m *model) impactAt(pos int, refBase, alt byte) Impact {
	key := pos<<8 | int(alt)
	if impact, ok := m.impacts[key]; ok {
		return impact
	}
	impact := m.ref.CodingImpact(m.contigs()[m.contig].Name, pos, string(refBase), string(alt))
	m.impacts[key] = impact
	return impact
}

// showReads is false once a read would be thinner than one column, where a
// pileup says nothing and the depth sketch says everything.
// readZoomLimit is the widest zoom that still draws individual reads. Past it
// a read is only a few columns wide and the depth track says the same thing
// for a fraction of the work.
const readZoomLimit = 400

// baseDetailSpan is how many reference bases may be on screen before the view
// stops keeping a base for every position. Beyond it no single base is drawn,
// so only the differences are worth keeping.
const baseDetailSpan = 10000

func (m *model) showReads() bool {
	return m.view.bpp <= min(max(m.readLength, 1), readZoomLimit)
}

func (m *model) measureReads() {
	total := 0
	for _, a := range m.reads {
		total += a.End - a.Start
	}
	if len(m.reads) > 0 {
		m.readLength = total / len(m.reads)
	}
}

func clamp(v, low, high int) int { return max(low, min(v, high)) }
