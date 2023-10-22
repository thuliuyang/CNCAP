/*
Copyright © 2017 the InMAP authors.
This file is part of InMAP.

InMAP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

InMAP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with InMAP.  If not, see <http://www.gnu.org/licenses/>.
*/

package inmap

import (
	"fmt"
	"math"
)

// Mech is an example type fulfils the github.com/spatialmodel/inmap/Mechanism
// interface.
type Mech struct{}

// Len returns the number of chemical species in this mechanism (9).
func (m Mech) Len() int {
	return 9
}

// emisConv lists the accepted names for emissions species, the array
// indices they correspond to, and the
// factors needed to convert [μg/s] of emitted species to [μg/s] of
// model species.
var emisConv = map[string]struct {
	i    int
	conv float64
}{
	"VOC":   {i: igOrg, conv: 1},
	"NOx":   {i: igNO, conv: NOxToN},
	"NH3":   {i: igNH, conv: NH3ToN},
	"SOx":   {i: igS, conv: SOxToS},
	"PM2_5": {i: iPM2_5, conv: 1},
}

// AddEmisFlux adds emissions flux to Cell c based on the given
// pollutant name and amount in units of μg/s. The units of
// the resulting flux are μg/m3/s.
func (m Mech) AddEmisFlux(c *Cell, name string, val float64) error {
	fluxScale := 1. / c.Dx / c.Dy / c.Dz // μg/s /m/m/m = μg/m3/s
	conv, ok := emisConv[name]
	if !ok {
		return fmt.Errorf("simplechem: '%s' is not a valid emissions species; valid options are VOC, NOx, NH3, SOx, and PM2_5", name)
	}
	c.EmisFlux[conv.i] += val * conv.conv * fluxScale
	return nil
}

// DryDep returns a dry deposition function of the type indicated by
// name that is compatible with this chemical mechanism.
// No options are available in this example.
func (m Mech) DryDep(name string) (CellManipulator, error) {
	panic("DryDep not implemented")
}

// WetDep returns a dry deposition function of the type indicated by
// name that is compatible with this chemical mechanism.
// No options are available in this example.
func (m Mech) WetDep(name string) (CellManipulator, error) {
	panic("WetDep not implemented")
}

// Species returns the names of the emission and concentration pollutant
// species that are used by this chemical mechanism.
func (m Mech) Species() []string {
	return []string{
		"VOCEmissions",
		"NOxEmissions",
		"NH3Emissions",
		"SOxEmissions",
		"PM25Emissions",
		"TotalPM25",
		"VOC",
		"SOA",
		"PrimaryPM25",
		"NH3",
		"pNH4",
		"SOx",
		"pSO4",
		"NOx",
		"pNO3",
	}
}

var emisLabels = map[string]int{"VOCEmissions": igOrg,
	"NOxEmissions":  igNO,
	"NH3Emissions":  igNH,
	"SOxEmissions":  igS,
	"PM25Emissions": iPM2_5,
}

// polLabels are labels and conversions for InMAP pollutants.
var polLabels = map[string]struct {
	index      []int     // index in concentration array
	conversion []float64 // conversion from N to NH4, S to SO4, etc...
}{
	"TotalPM25": {[]int{iPM2_5, ipOrg, ipNH, ipS, ipNO},
		[]float64{1, 1, NtoNH4, StoSO4, NtoNO3}},
	"VOC":         {[]int{igOrg}, []float64{1.}},
	"SOA":         {[]int{ipOrg}, []float64{1.}},
	"PrimaryPM25": {[]int{iPM2_5}, []float64{1.}},
	"NH3":         {[]int{igNH}, []float64{1. / NH3ToN}},
	"pNH4":        {[]int{ipNH}, []float64{NtoNH4}},
	"SOx":         {[]int{igS}, []float64{1. / SOxToS}},
	"pSO4":        {[]int{ipS}, []float64{StoSO4}},
	"NOx":         {[]int{igNO}, []float64{1. / NOxToN}},
	"pNO3":        {[]int{ipNO}, []float64{NtoNO3}},
}

// Value returns the concentration or emissions value of
// the given variable in the given Cell. It returns an
// error if given an invalid variable name.
func (m Mech) Value(c *Cell, variable string) (float64, error) {
	i, ok := emisLabels[variable]
	if ok {
		return c.EmisFlux[i], nil
	}
	conv, ok := polLabels[variable]
	if !ok {
		return math.NaN(), fmt.Errorf("simplechem: invalid variable name %s; valid names are %v", variable, m.Species())
	}
	var val float64
	for ii, i := range conv.index {
		val += c.Cf[i] * conv.conversion[ii]
	}
	return val, nil
}

// Units returns the units of the given variable, or an
// error if the variable name is invalid.
func (m Mech) Units(variable string) (string, error) {
	if _, ok := emisLabels[variable]; ok {
		return "μg/m³/s", nil
	}
	if _, ok := polLabels[variable]; !ok {
		return "", fmt.Errorf("simplechem: invalid variable name %s; valid names are %v", variable, m.Species())
	}
	return "μg/m³", nil
}

// Chemistry returns a function that calculates the secondary formation of PM2.5.
// It explicitly calculates formation of particulate sulfate
// from gaseous and aqueous SO2.
// It partitions organic matter ("gOrg" and "pOrg"), the
// nitrogen in nitrate ("gNO and pNO"), and the nitrogen in ammonia ("gNH" and
// "pNH) between gaseous and particulate phase
// based on the spatially explicit partioning present in the baseline data.
// The function arguments represent the array indices of each chemical species.
func (m Mech) Chemistry() CellManipulator {
	return func(c *Cell, Δt float64) {
		// All SO4 forms particles, so sulfur particle formation is limited by the
		// SO2 -> SO4 reaction.
		ΔS := c.SO2oxidation * c.Cf[igS] * Δt
		c.Cf[ipS] += ΔS
		c.Cf[igS] -= ΔS
		// NH3 / pNH4 partitioning
		totalNH := c.Cf[igNH] + c.Cf[ipNH]
		c.Cf[ipNH] = totalNH * c.NHPartitioning
		c.Cf[igNH] = totalNH * (1 - c.NHPartitioning)

		// NOx / pN0 partitioning
		totalNO := c.Cf[igNO] + c.Cf[ipNO]
		c.Cf[ipNO] = totalNO * c.NOPartitioning
		c.Cf[igNO] = totalNO * (1 - c.NOPartitioning)

		// VOC/SOA partitioning
		totalOrg := c.Cf[igOrg] + c.Cf[ipOrg]
		c.Cf[ipOrg] = totalOrg * c.AOrgPartitioning
		c.Cf[igOrg] = totalOrg * (1 - c.AOrgPartitioning)
	}
}
