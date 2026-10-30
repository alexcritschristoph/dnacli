package main

import (
	"fmt"
	"os"
	"testing"
)

// TestMain keeps the tests off the clipboard of whoever is running them, and
// skips the suite when the fixtures in testdata/ are not here. The datasets add
// up to 429 MB, so they are not in the repository: a clone without them runs
// nothing, rather than failing on files it has no way to read.
func TestMain(m *testing.M) {
	clipboardTools = nil
	os.Unsetenv("TMUX")
	os.Unsetenv("SSH_TTY")
	if _, err := os.Stat("testdata"); err != nil {
		fmt.Println("Skipping every test: the fixtures in testdata/ are not here.")
		os.Exit(0)
	}
	os.Exit(m.Run())
}
