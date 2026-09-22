package pde

import (
	"bytes"
	"math"
	"testing"
)

func mustMarshal(t *testing.T, s *Snapshot) []byte {
	t.Helper()
	b, err := s.Marshal()
	if err != nil {
		t.Fatal(err)
	}
	return b
}

func nan() float64 { return math.NaN() }

// injectToken replaces the first numeric token inside the JSON "data" array
// with token (e.g. Infinity), simulating a corrupted checkpoint file.
func injectToken(raw []byte, token string) []byte {
	key := []byte(`"data": [`)
	start := bytes.Index(raw, key)
	if start < 0 {
		panic("data key not found")
	}
	start += len(key)
	for start < len(raw) && (raw[start] == ' ' || raw[start] == '\n' || raw[start] == '\t') {
		start++
	}
	end := start
	for end < len(raw) && raw[end] != ',' && raw[end] != ']' {
		end++
	}
	out := bytes.Clone(raw[:start])
	out = append(out, []byte(token)...)
	out = append(out, raw[end:]...)
	return out
}
