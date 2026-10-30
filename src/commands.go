package main

import (
	"strconv"
	"strings"
	"unicode"

	tea "github.com/charmbracelet/bubbletea"
)

// A command is something the user can type into the bar along the top.
type command struct {
	name   string
	help   string
	needs  string // what to type after the name, when it wants an argument
	hidden bool   // a spelling that works but is not listed
	run    func(m *model, arg string)
}

var commands []command

// init fills the command list. It is not a var initializer because the export
// command draws a frame, and drawing a frame reads this list for suggestions.
func init() {
	commands = []command{
		{name: "contigs", help: "show contigs in BAM", run: func(m *model, _ string) { m.openPanel("contigs") }},
		{name: "colors", help: "change read highlighting", run: func(m *model, _ string) { m.openPanel("colors") }},
		{name: "colours", hidden: true, run: func(m *model, _ string) { m.openPanel("colors") }},
		{name: "filters", help: "filter reads in file", run: (*model).runFilter},
		{name: "variants", help: "show variants in VCF", run: func(m *model, _ string) { m.openPanel("variants") }},
		{name: "go", help: "go to a position on this contig", needs: "position", run: (*model).parseGoto},
		{name: "gene", help: "go to a gene by name", needs: "name", run: (*model).goToGene},
		{name: "open", help: "open a BAM, GenBank or VCF file", needs: "path", run: (*model).runOpen},
		{name: "export", help: "save this view as an SVG", needs: "[file] [light|dark]", run: (*model).runExport},
	}
}

// suggestion is one line of the list under the command bar.
type suggestion struct {
	label   string
	help    string
	replace string // what the bar becomes when this is chosen; a trailing space means more to type
}

// suggestions are what the text in the bar could become.
func (m *model) suggestions() []suggestion {
	name, arg, typedArg := strings.Cut(strings.TrimLeft(m.command, " "), " ")
	if typedArg {
		switch strings.ToLower(name) {
		case "gene":
			return m.geneSuggestions(arg)
		case "filters", "filter":
			return filterSuggestions(arg)
		case "go":
			c := m.contigs()[m.contig]
			middle := c.Length / 2
			return []suggestion{{comma(middle), "position", "go "},
				{c.Name + ":" + comma(middle), "contig: position", "go "},
				{comma(max(1, middle-100)) + "-" + comma(min(c.Length, middle+200)), "a region", "go "}}
		case "open":
			if arg == "" {
				return []suggestion{
					{"reads.bam", "another BAM file in a new tab", "open "},
					{"genes.gbk", "a GenBank reference", "open "},
					{"calls.vcf", "a VCF of variants", "open "},
				}
			}
		case "export":
			return []suggestion{
				{"figure.svg light", "save view as SVG (light)", "export figure.svg light"},
				{"figure.svg dark", "save view as SVG (dark)", "export figure.svg dark"},
			}
		}
		return nil
	}
	var out []suggestion
	for _, c := range commands {
		if c.hidden || !strings.HasPrefix(c.name, strings.ToLower(name)) {
			continue
		}
		replace := c.name
		if c.needs != "" {
			replace += " "
		}
		out = append(out, suggestion{c.name, c.help, replace})
	}
	return out
}

func filterSuggestions(typed string) []suggestion {
	all := []suggestion{
		{"mapq 30", "minimum mapq", "filters mapq "},
		{"quality 20", "minimum average base quality", "filters quality "},
		{"mismatches 2", "maximum mismatches from reference", "filters mismatches "},
		{"length 500", "minimum read length", "filters length "},
		{"secondary", "show or hide secondary alignments", "filters secondary"},
		{"duplicates", "show or hide duplicates", "filters duplicates"},
		{"pairs", "show or hide reads not in a proper pair", "filters pairs"},
		{"clips", "show or hide soft-clipped bases", "filters clips"},
		{"reset", "reset all filters", "filters reset"},
	}
	var out []suggestion
	for _, s := range all {
		if strings.HasPrefix(s.label, strings.ToLower(typed)) {
			out = append(out, s)
		}
	}
	return out
}

func (m *model) geneSuggestions(typed string) []suggestion {
	var out []suggestion
	for _, g := range m.genes(typed) {
		where := g.contig + ":" + comma(g.start+1)
		if len(out) < 8 {
			out = append(out, suggestion{g.name, where, "gene " + g.name})
		}
	}
	return out
}

// gene is one annotation the user can jump to by name.
type gene struct {
	contig, name string
	start, end   int
}

// genes lists the coding sequences whose name starts with what was typed, the
// contig on screen first.
func (m *model) genes(prefix string) []gene {
	if m.ref == nil {
		return nil
	}
	var out []gene
	seen := map[string]bool{}
	for i := range m.contigs() {
		c := m.contigs()[(m.contig+i)%len(m.contigs())]
		for _, f := range m.ref.Features[c.Name] {
			name := geneName(f)
			if f.Kind != "CDS" || seen[c.Name+name] || !strings.HasPrefix(strings.ToLower(name), strings.ToLower(prefix)) {
				continue
			}
			seen[c.Name+name] = true
			out = append(out, gene{c.Name, name, f.Start(), f.End()})
		}
	}
	return out
}

// commandRune is true for the characters commands are made of: letters,
// digits, and the punctuation of positions, gene names and file paths. Any
// other key belongs to the view rather than the bar, as does anything a
// terminal sends that is not plain text.
func commandRune(r rune) bool {
	return r < unicode.MaxASCII && (unicode.IsLetter(r) || unicode.IsDigit(r) || strings.ContainsRune(" ,:._-/~", r))
}

// startsCommand is true when text is the start of a command, or a whole command
// with its argument, as a burst of quick typing can be.
// Only text like that opens the bar, because a terminal kept busy by the mouse
// sends stray characters of its own between the reports it means to send.
func startsCommand(text string) bool {
	for _, c := range commands {
		if strings.HasPrefix(c.name, text) || strings.HasPrefix(text, c.name) {
			return true
		}
	}
	return false
}

// editCommand adds a keystroke to the command bar and says whether the bar took
// it. Typed text arrives a burst at a time when someone types quickly, or
// pastes.
func (m *model) editCommand(msg tea.KeyMsg) bool {
	switch msg.Type {
	case tea.KeyRunes:
		text := string(msg.Runes)
		if msg.Alt || (m.command == "" && !startsCommand(text)) {
			return false
		}
		probe := m.command
		for _, r := range msg.Runes {
			if !commandRune(r) {
				return false
			}
			probe += string(r)
		}
		m.command, m.choice = probe, 0
		return true
	case tea.KeySpace:
		if m.command == "" {
			return false // a space does not start a command
		}
		m.command, m.choice = m.command+" ", 0
		return true
	}
	switch msg.String() {
	case "esc":
		m.command, m.choice = "", 0
	case "enter":
		m.enterCommand()
	case "tab":
		m.completeCommand()
	case "up":
		m.choice = max(0, m.choice-1)
	case "down":
		m.choice = min(m.choice+1, max(0, len(m.suggestions())-1))
	case "backspace":
		if r := []rune(m.command); len(r) > 0 {
			m.command, m.choice = string(r[:len(r)-1]), 0
		}
	default:
		return false
	}
	return true
}

// enterCommand runs the line picked out of the list. A command that wants
// more is filled in instead, and the list turns into what could follow it.
func (m *model) enterCommand() {
	typed := m.command
	m.completeCommand()
	if m.command != typed && strings.HasSuffix(m.command, " ") {
		return
	}
	m.runCommand()
}

// completeCommand fills the bar in with the highlighted line, unless more has
// been typed than the line would leave.
func (m *model) completeCommand() {
	list := m.suggestions()
	if len(list) == 0 {
		return
	}
	if pick := list[clamp(m.choice, 0, len(list)-1)].replace; strings.HasPrefix(strings.ToLower(pick), strings.ToLower(m.command)) {
		m.command, m.choice = pick, 0
	}
}

// runCommand does what the bar says, and empties it.
func (m *model) runCommand() {
	typed := strings.TrimSpace(m.command)
	m.command, m.choice = "", 0
	name, arg, _ := strings.Cut(typed, " ")
	if name == "" {
		return
	}
	if _, err := strconv.Atoi(strings.ReplaceAll(name, ",", "")); err == nil {
		m.parseGoto(typed) // a bare number is a position
		return
	}
	for _, c := range commands {
		if strings.HasPrefix(c.name, strings.ToLower(name)) {
			c.run(m, strings.TrimSpace(arg))
			return
		}
	}
	m.notify("There is no %q command. Type to see what there is.", name)
}

// runFilter changes one filter and shows the panel, so that the change is
// visible next to the others.
func (m *model) runFilter(arg string) {
	word, value, _ := strings.Cut(arg, " ")
	switch strings.ToLower(word) {
	case "":
	case "mapq":
		m.setFilterNumber(&m.filter.MinMapQ, value, 0, 60)
	case "quality":
		m.setFilterNumber(&m.filter.MinMeanQual, value, 0, 40)
	case "mismatches":
		if strings.EqualFold(value, "any") {
			m.filter.MaxMismatches = -1
		} else {
			m.setFilterNumber(&m.filter.MaxMismatches, value, 0, 50)
		}
	case "length":
		m.setFilterNumber(&m.filter.MinLength, value, 0, 1000000)
	case "secondary":
		m.filter.HideSecondary = !m.filter.HideSecondary
	case "duplicates":
		m.filter.HideDuplicates = !m.filter.HideDuplicates
	case "pairs":
		m.filter.HideUnpaired = !m.filter.HideUnpaired
	case "clips":
		m.filter.SoftClips = !m.filter.SoftClips
	case "reset":
		m.filter = defaultFilter()
	default:
		m.notify("There is no %q filter.", word)
	}
	m.openPanel("filters")
}

func (m *model) setFilterNumber(field *int, value string, low, high int) {
	n, err := strconv.Atoi(value)
	if err != nil {
		m.notify("%q is not a number.", value)
		return
	}
	*field = clamp(n, low, high)
}

// goToGene frames the gene with that name, wherever it is.
func (m *model) goToGene(name string) {
	if name == "" {
		m.notify("Type part of a gene name after gene.")
		return
	}
	found := m.genes(name)
	if len(found) == 0 {
		m.notify("No gene called %q.", name)
		return
	}
	g := found[0]
	for i, c := range m.contigs() {
		if c.Name == g.contig && i != m.contig {
			m.selectContig(i)
		}
	}
	m.showRange(g.start, g.end)
}
