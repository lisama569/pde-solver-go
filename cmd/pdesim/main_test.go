package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/lisama569/pde-solver-go/pde"
)

// TestCLISplitResume drives the real run() entry point: a partial run writes a
// checkpoint, a resume invocation (different workers, extended total-time)
// finishes, and the result equals an uninterrupted serial run.
func TestCLISplitResume(t *testing.T) {
	dir := t.TempDir()
	ck := filepath.Join(dir, "ck.json")
	base := []string{
		"-nx=18", "-ny=14", "-dx=0.05", "-dy=0.05",
		"-alpha=0.15", "-dt=0.002", "-x0=d=1.5", "-x1=n=0",
		"-y0=d=0", "-y1=n=0.2", "-source-scale=2.5", "-source-cx=0.4",
	}

	out := capture(func() {
		if err := run(append(append([]string{"run"}, base...),
			"-total-time=0.01", "-workers=2", "-checkpoint="+ck)); err != nil {
			t.Fatal(err)
		}
	})
	partialSum := parseSum(t, out)

	out = capture(func() {
		if err := run(append(append([]string{"resume"}, base...),
			"-total-time=0.02", "-workers=5", "-checkpoint="+ck,
			"-check-invariant")); err != nil {
			t.Fatal(err)
		}
	})
	resumedSum := parseSum(t, out)

	out = capture(func() {
		if err := run(append(append([]string{"run"}, base...),
			"-total-time=0.02", "-workers=1")); err != nil {
			t.Fatal(err)
		}
	})
	continuousSum := parseSum(t, out)

	t.Logf("partial sum=%s resumed sum=%s continuous sum=%s", partialSum, resumedSum, continuousSum)
	if resumedSum != continuousSum {
		t.Fatalf("resumed (%s) != continuous (%s)", resumedSum, continuousSum)
	}
	if _, err := os.Stat(ck); err != nil {
		t.Fatalf("checkpoint missing: %v", err)
	}
}

// TestCLIStabilityFailureTyped confirms the entry point surfaces the typed
// stability failure (distinct from config errors) and exits via run error.
func TestCLIStabilityFailureTyped(t *testing.T) {
	err := run([]string{
		"run", "-nx=10", "-ny=10", "-dx=0.1", "-dy=0.1", "-alpha=1",
		"-dt=0.1", "-total-time=0.1",
	})
	f := pde.AsFailure(err)
	if f == nil || f.Kind != pde.KindStability {
		t.Fatalf("want KindStability from CLI, got %v", err)
	}
	t.Logf("CLI stability failure: %s", f.Msg)
	if !strings.Contains(f.Msg, "dt must be <=") {
		t.Fatalf("stability message should state dt bound: %s", f.Msg)
	}
}

func capture(fn func()) string {
	old := os.Stdout
	r, w, _ := os.Pipe()
	os.Stdout = w
	fn()
	w.Close()
	os.Stdout = old
	b, _ := readAll(r)
	return string(b)
}

func readAll(f *os.File) ([]byte, error) {
	buf := make([]byte, 0, 4096)
	tmp := make([]byte, 1024)
	for {
		n, err := f.Read(tmp)
		buf = append(buf, tmp[:n]...)
		if err != nil {
			return buf, nil
		}
	}
}

func parseSum(t *testing.T, out string) string {
	t.Helper()
	for _, line := range strings.Split(out, "\n") {
		if strings.HasPrefix(line, "steps=") {
			fields := strings.Fields(line)
			return fields[len(fields)-1]
		}
	}
	t.Fatalf("no steps= line in output: %q", out)
	return ""
}
