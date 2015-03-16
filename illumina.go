// Copyright ©2013 The bíogo.illumina Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package illumina provides support for handling Illumina read metadata.
package illumina

import (
	"errors"
	"fmt"
	"strconv"
	"strings"

	"github.com/biogo/biogo/alphabet"
)

var (
	ErrBadIdentifer = errors.New("illumina: unable to parse identifier")
	ErrBadTag       = errors.New("illumina: illegal multiplex tag")
)

type Interface interface {
	Name() string
	Description() string
}

type Type int

func (t Type) String() string {
	switch t {
	case Undefined:
		return "undefined"
	case PreCasava:
		return "pre-casava"
	case Casava:
		return "casava"
	}
	return "invalid"
}

const (
	Undefined Type = iota
	PreCasava
	Casava
)

// A Coordinate represents a cluster location in an Illumina flow-cell lane.
type Coordinate struct {
	X, Y int
}

// A Multiplex represents multiplexing tag information.
type Multiplex struct {
	Index int8   // Index is -1 if not valid.
	Tag   string // Tag is empty if not valid.
}

// A Metadata represents Illumina read metadata.
type Metadata struct {
	Type        Type
	Instrument  string     // Unique instrument name.
	Run         int        // Run id, -1 if not valid.
	FlowCell    string     // Flowcell id.
	Lane        int8       // Flowcell lane.
	Tile        int        // Tile number within the flowcell lane.
	Coordinate  Coordinate // Coordinate of the cluster within the tile.
	Mate        int8       // Member of a pair, 1 or 2 for paired reads.
	BadRead     bool       // Read failed filter.
	ControlBits int        // 0 when none of the control bits are on, otherwise it is an even number, -1 if not valid.
	Multiplex   Multiplex  // Multiplexing information.
}

// Parse parses the name and description fields of an Illumina identifier represented by r.
// An error is returned if the input does not conform to the pre-Casava or Casava identifier
// formats.
func Parse(r Interface) (Metadata, error) {
	name := r.Name()
	if strings.Index(name, "#") >= 0 {
		return preCasava(name)
	}
	desc := r.Description()
	return casava(name, desc)
}

func mustAtoi(s string) int {
	if len(s) == 0 {
		return -1
	}
	i, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return i
}

func atob(s string) (int8, error) {
	if len(s) == 0 {
		return -1, nil
	}
	i, err := strconv.Atoi(s)
	return int8(i), err
}

func tagOk(tag string) bool {
	for _, r := range tag {
		if !alphabet.DNA.IsValid(alphabet.Letter(r)) {
			return false
		}
	}
	return true
}

// @HWUSI-EAS100R:6:73:941:1973#0/1
//
//  HWUSI-EAS100R 	the unique instrument name
//  6 				flowcell lane
//  73 				tile number within the flowcell lane
//  941 			'x'-coordinate of the cluster within the tile
//  1973		 	'y'-coordinate of the cluster within the tile
//  #0 				index number for a multiplexed sample (0 for no indexing)
//  /1 				the member of a pair, /1 or /2 (paired-end or mate-pair reads only)
func preCasavaSep(r rune) bool { return r == ':' || r == '#' || r == '/' }
func preCasava(name string) (m Metadata, err error) {
	f := strings.FieldsFunc(name, preCasavaSep)
	if len(f) < 6 {
		return Metadata{}, ErrBadIdentifer
	}
	defer func() {
		r := recover()
		if e, ok := r.(error); ok {
			err = e
			m.Type = Undefined
		}
	}()
	m = Metadata{
		Type:       PreCasava,
		Instrument: f[0],
		Run:        -1,
		Lane:       int8(mustAtoi(f[1])),
		Tile:       mustAtoi(f[2]),
		Coordinate: Coordinate{mustAtoi(f[3]), mustAtoi(f[4])},
	}
	m.Multiplex.Index, err = atob(f[5])
	if err != nil {
		if tagOk(f[5]) {
			err = nil
			m.Multiplex.Tag = f[5]
			m.Multiplex.Index = -1
		} else {
			panic(fmt.Errorf("illumina: illegal multiplex tag: %s", err))
		}
	}
	if len(f) == 7 {
		m.Mate = int8(mustAtoi(f[6]))
	}
	return m, err
}

// @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
//
// EAS139 	the unique instrument name
// 136 		the run id
// FC706VJ 	the flowcell id
// 2 		flowcell lane
// 2104 	tile number within the flowcell lane
// 15343 	'x'-coordinate of the cluster within the tile
// 197393 	'y'-coordinate of the cluster within the tile
// 1 		the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
// Y 		Y if the read fails filter (read is bad), N otherwise
// 18	 	0 when none of the control bits are on, otherwise it is an even number
// ATCACG 	index sequence
func casavaSep(r rune) bool { return r == ':' }
func casava(name, desc string) (m Metadata, err error) {
	nf := strings.FieldsFunc(name, casavaSep)
	df := strings.FieldsFunc(desc, casavaSep)
	if !(len(nf) == 7 && (len(df) == 4 || desc == "")) {
		return Metadata{}, ErrBadIdentifer
	}
	if len(df) == 4 && !tagOk(df[3]) {
		return Metadata{}, ErrBadTag
	}
	defer func() {
		r := recover()
		if e, ok := r.(error); ok {
			err = e
			m.Type = Undefined
		}
	}()
	m = Metadata{
		Type:        Casava,
		Instrument:  nf[0],
		Run:         mustAtoi(nf[1]),
		FlowCell:    nf[2],
		Lane:        int8(mustAtoi(nf[3])),
		Tile:        mustAtoi(nf[4]),
		Coordinate:  Coordinate{mustAtoi(nf[5]), mustAtoi(nf[6])},
		ControlBits: -1,
		Multiplex:   Multiplex{Index: -1},
	}
	if desc != "" {
		m.Mate = int8(mustAtoi(df[0]))
		m.BadRead = df[1] == "Y" || df[1] == "y"
		m.ControlBits = mustAtoi(df[2])
		m.Multiplex.Tag = df[3]
	}

	return m, nil
}
