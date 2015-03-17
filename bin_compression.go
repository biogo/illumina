// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package illumina

import (
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/quality"
)

type Scheme *[256]alphabet.Qphred

// DefaultCompression uses the quality compression scheme described in the Illumina
// white paper :
//   Old Quality Score   New Quality score
//         0-1                   0
//         2–9                   6
//        10–19                 15
//        20–24                 22
//        25–29                 27
//        30–34                 33
//        35–39                 37
//         ≥ 40                 40
var DefaultCompression Scheme = defaultCompression

var defaultCompression = func() Scheme {
	var cs [256]alphabet.Qphred
	for i := 2; i < 256; i++ {
		var q alphabet.Qphred
		switch {
		case 2 <= i && i < 10:
			q = 6
		case 10 <= i && i < 20:
			q = 15
		case 20 <= i && i < 25:
			q = 22
		case 25 <= i && i < 30:
			q = 27
		case 30 <= i && i < 35:
			q = 33
		case 35 <= i && i < 40:
			q = 37
		case i >= 40:
			q = 40
		}
		cs[i] = q
	}
	return &cs
}()

// BinCompress lossily compresses the qualities of a seq.Scorer according to the
// the provided compression Scheme, c. If c is nil, the default Scheme is used. The
// approach used by BinCompress is described in the  Illumina whitepaper at
// http://www.illumina.com/Documents/products/whitepapers/whitepaper_datacompression.pdf
// with the exception that N is presumed from based quality.
// Solexa qualities are handled by translating via a calculated E value.
func BinCompress(s seq.Scorer, c Scheme) error {
	if c == nil {
		c = DefaultCompression
	}
	if sl, ok := s.(seq.Slicer); ok {
		switch d := sl.Slice().(type) {
		case alphabet.QLetters:
			for i, ql := range d {
				d[i].Q = c[ql.Q]
			}
		case quality.Qphreds:
			for i, q := range d {
				d[i] = c[q]
			}
		default:
			return slowCompression(s, c)
		}
		return nil
	}
	return slowCompression(s, c)
}

func slowCompression(s seq.Scorer, c Scheme) error {
	for i := s.Start(); i < s.End(); i++ {
		err := s.SetE(i, (*c)[alphabet.Ephred(s.EAt(i))].ProbE())
		if err != nil {
			return err
		}
	}
	return nil
}
