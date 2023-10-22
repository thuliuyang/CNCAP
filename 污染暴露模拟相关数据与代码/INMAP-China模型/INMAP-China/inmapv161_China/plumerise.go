/*
Copyright © 2013 the InMAP authors.
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
//	"github.com/ctessum/atmos/plumerise"
	"fmt"
	"errors"
	"math"
	"sort"
)
// IsPlumeIn calculates whether the plume rise from an emission is at the height
// of c when given stack information
// (see github.com/ctessum/atmos/plumerise for required units).
// The return values are whether the plume rise ends within the current cell,
// the height of the plume rise in meters, and whether there was an error.
func (c *Cell) IsPlumeIn(stackHeight, stackDiam, stackTemp, stackVel float64) (bool, float64, error) {

	// Find the cells in the vertical column below c.
	var cellStack []*Cell
	cc := c
	for {
		cellStack = append(cellStack, cc)
		if (*cc.groundLevel)[0].Cell == cc {
			break
		}
		cc = (*cc.below)[0].Cell
	}
	// reverse the order of the stack so it starts at ground level.
	for left, right := 0, len(cellStack)-1; left < right; left, right = left+1, right-1 {
		cellStack[left], cellStack[right] = cellStack[right], cellStack[left]
	}

	layerHeights := make([]float64, len(cellStack)+1)
	temperature := make([]float64, len(cellStack))
	windSpeed := make([]float64, len(cellStack))
	windSpeedInverse := make([]float64, len(cellStack))
	windSpeedMinusThird := make([]float64, len(cellStack))
	windSpeedMinusOnePointFour := make([]float64, len(cellStack))
	sClass := make([]float64, len(cellStack))
	s1 := make([]float64, len(cellStack))

	for i, cell := range cellStack {
		layerHeights[i+1] = layerHeights[i] + cell.Dz
		temperature[i] = cell.Temperature
		windSpeed[i] = cell.WindSpeed
		windSpeedInverse[i] = cell.WindSpeedInverse
		windSpeedMinusThird[i] = cell.WindSpeedMinusThird
		windSpeedMinusOnePointFour[i] = cell.WindSpeedMinusOnePointFour
		sClass[i] = cell.SClass
		s1[i] = cell.S1
	}

	plumeIndex, plumeHeight, err := ASMEPrecomputed(stackHeight, stackDiam,
		stackTemp, stackVel, layerHeights, temperature, windSpeed,
		sClass, s1, windSpeedMinusOnePointFour, windSpeedMinusThird,
		windSpeedInverse)

	if err != nil {
		if err == ErrAboveModelTop {
			// If the plume is above the top of our stack, return true if c is
			// in the top model layer (because we want to put the plume in the
			// top layer even if it should technically go above it),
			//  otherwise return false.
			if (*c.above)[0].boundary {
				return true, plumeHeight, nil
			}
			return false, plumeHeight, nil
		}
		return false, plumeHeight, err
	}

	// if the index of the plume is at the end of the cell stack,
	// that means that the plume should go in this cell.
	if plumeIndex == c.Layer {
		return true, plumeHeight, nil
	}
	return false, plumeHeight, nil
}

func ASME(stackHeight, stackDiam, stackTemp,
	stackVel float64, layerHeights, temperature, windSpeed,
	sClass, s1 []float64) (plumeLayer int, plumeHeight float64, err error) {

	stackLayer, err := findLayer(layerHeights, stackHeight)
	if err != nil {
		return stackLayer, stackHeight, err
	}
	deltaH, err := calcDeltaH(stackLayer, temperature, windSpeed, sClass, s1,
		stackHeight, stackTemp, stackVel, stackDiam)
	if err != nil {
		return
	}

	plumeHeight = stackHeight + deltaH
	plumeLayer, err = findLayer(layerHeights, plumeHeight)
	return
}

// ASMEPrecomputed is the same as ASME except it takes
// precomputed (averaged) meteorological parameters:
// the inverse of the stability parameter (s1Inverse [1/unknown units],
// unstaggered grid), windSpeedMinusOnePointFour [(m/s)^(-1.4)] (unstaggered grid),
// windSpeedMinusThird [(m/s)^(-1/3)] (unstaggered grid),
// and windSpeedInverse [(m/s)^(-1)] (unstaggered grid),
// Uses the plume rise calculation: ASME (1973), as described in Sienfeld and Pandis,
// ``Atmospheric Chemistry and Physics - From Air Pollution to Climate Change
func ASMEPrecomputed(stackHeight, stackDiam, stackTemp,
	stackVel float64, layerHeights, temperature, windSpeed,
	sClass, s1, windSpeedMinusOnePointFour, windSpeedMinusThird,
	windSpeedInverse []float64) (plumeLayer int, plumeHeight float64, err error) {

	stackLayer, err := findLayer(layerHeights, stackHeight)
	if err != nil {
		return stackLayer, stackHeight, err
	}
	deltaH, err := calcDeltaHPrecomputed(stackLayer, temperature, windSpeed, sClass,
		s1, stackHeight, stackTemp, stackVel, stackDiam,
		windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse)
	if err != nil {
		return
	}

	plumeHeight = stackHeight + deltaH
//	fmt.Println("PH")
//	fmt.Println(plumeHeight)
	plumeLayer, err = findLayer(layerHeights, plumeHeight)
	return
}

// Find K level of stack or plume
func findLayer(layerHeights []float64, stackHeight float64) (int, error) {
	stackLayer := sort.SearchFloat64s(layerHeights, stackHeight)
	if stackLayer == len(layerHeights) {
		stackLayer -= 2
		return stackLayer, ErrAboveModelTop
	}
	if stackLayer != 0 {
		stackLayer--
	}
	return stackLayer, nil
}

// calcDeltaH calculates plume rise (ASME, 1973).
func calcDeltaH(stackLayer int, temperature, windSpeed, sClass, s1 []float64,
	stackHeight, stackTemp, stackVel, stackDiam float64) (float64, error) {
	deltaH := 0. // Plume rise, (m).

	airTemp := temperature[stackLayer]
	windSpd := windSpeed[stackLayer]

	if (stackTemp-airTemp) < 50. &&
		stackVel > windSpd && stackVel > 10. {

		// Plume is dominated by momentum forces
		deltaH = stackDiam * math.Pow(stackVel, 1.4) / math.Pow(windSpd, 1.4)

	} else { // Plume is dominated by buoyancy forces

		// Bouyancy flux, m4/s3
		F := g * (stackTemp - airTemp) / stackTemp * stackVel *
			math.Pow(stackDiam/2, 2)

		if sClass[stackLayer] > 0.5 { // stable conditions

			deltaH = 29. * math.Pow(
				F/s1[stackLayer], 0.333333333) /
				math.Pow(windSpd, 0.333333333)

		} else { // unstable conditions

			deltaH = 7.4 * math.Pow(F*math.Pow(stackHeight, 2.),
				0.333333333) / windSpd
		}
	}
	if math.IsNaN(deltaH) {
		err := fmt.Errorf("plume height == NaN\n"+
			"deltaH: %v, stackDiam: %v,\n"+
			"stackVel: %v, windSpd: %v, stackTemp: %v,\n"+
			"airTemp: %v, stackHeight: %v\n",
			deltaH, stackDiam, stackVel,
			windSpd, stackTemp, airTemp, stackHeight)
		return deltaH, err
	}
	return deltaH, nil
}

// calcDeltaHPrecomputed calculates plume rise, the same as calcDeltaH,
// (ASME, 1973), except that it uses precomputed meteorological parameters.
func calcDeltaHPrecomputed(stackLayer int, temperature, windSpeed, sClass,
	s1 []float64,
	stackHeight, stackTemp, stackVel, stackDiam float64,
	windSpeedMinusOnePointFour, windSpeedMinusThird,
	windSpeedInverse []float64) (float64, error) {

	deltaH := 0. // Plume rise, (m).

	airTemp := temperature[stackLayer]
	windSpd := windSpeed[stackLayer]

	if (stackTemp-airTemp) < 50. &&
		stackVel > windSpd && stackVel > 10. {

		// Plume is dominated by momentum forces
		deltaH = stackDiam * math.Pow(stackVel, 1.4) *
			windSpeedMinusOnePointFour[stackLayer]
        	      //fmt.Println("plume related1:")
	              //fmt.Println(deltaH,stackTemp,airTemp,stackHeight, stackDiam,stackVel,windSpeedMinusOnePointFour[stackLayer])
		if math.IsNaN(deltaH) {
			return deltaH, fmt.Errorf("plumerise: momentum-dominated deltaH is NaN. "+
				"stackDiam: %g, stackVel: %g, windSpeedMinusOnePointFour: %g",
				stackDiam, stackVel, windSpeedMinusOnePointFour[stackLayer])
		}

	} else { // Plume is dominated by buoyancy forces

		var tempDiff float64
		if stackTemp-airTemp == 0 {
			tempDiff = 0
		} else {
			tempDiff = 2 * (stackTemp - airTemp) / (stackTemp + airTemp)
		}

		// Bouyancy flux, m4/s3
		F := g * tempDiff * stackVel *
			math.Pow(stackDiam/2, 2)

		if sClass[stackLayer] > 0.5 && s1[stackLayer] != 0 && F > 0 { // stable conditions

			// Ideally, we would also use the inverse of S1,
			// but S1 is zero sometimes so that doesn't work well.
			deltaH = 29. * math.Pow(
				F/s1[stackLayer], 0.333333333) * windSpeedMinusThird[stackLayer]

			if math.IsNaN(deltaH) {
				return deltaH, fmt.Errorf("plumerise: stable bouyancy-dominated deltaH is NaN. "+
					"F: %g, s1: %g, windSpeedMinusThird: %g",
					F, s1[stackLayer], windSpeedMinusThird[stackLayer])
			}

		} else if F > 0. { // unstable conditions

			deltaH = 7.4 * math.Pow(F*math.Pow(stackHeight, 2.),
				0.333333333) * windSpeedInverse[stackLayer]

			if math.IsNaN(deltaH) {
				return deltaH, fmt.Errorf("plumerise: unstable bouyancy-dominated deltaH is NaN. "+
					"F: %g, stackHeight: %g, windSpeedInverse: %g",
					F, stackHeight, windSpeedInverse[stackLayer])
			}
		} else {
			// If F < 0, the unstable algorithm above will give an imaginary plume rise.
			deltaH = 0
		}
	//fmt.Println("plume related:")
	//fmt.Println(deltaH,F,stackTemp,airTemp,stackHeight, stackDiam,stackVel)
	}
	return deltaH, nil
}

// ErrAboveModelTop is returned when the plume is above the top
// model layer.
var ErrAboveModelTop = errors.New("plume rise > top of grid")

