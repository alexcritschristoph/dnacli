// Package lipgloss stands in for the real library, which bubbletea asks for
// the terminal's background colour at startup. That query blocks for five
// seconds on a terminal that does not answer it, and eats keys typed meanwhile.
// dnacli sets its own colours, so the answer is never used.
package lipgloss

// HasDarkBackground reports the dark background dnacli draws on.
func HasDarkBackground() bool { return true }
