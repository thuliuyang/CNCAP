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

package simpledrydep_test

import (
	"testing"

	"github.com/spatialmodel/inmap"
	"github.com/spatialmodel/inmap/science/chem/simplechem"
	"github.com/spatialmodel/inmap/science/drydep/simpledrydep"
)

// Indicies of individual pollutants in arrays.
const (
	igOrg int = iota
	ipOrg
	iPM2_5
	igNH
	ipNH
	igS
	ipS
	igNO
	ipNO
)

// simpleDryDepIndices provides array indices for use with package simpledrydep.
func simpleDryDepIndices() (simpledrydep.SOx, simpledrydep.NH3, simpledrydep.NOx, simpledrydep.VOC, simpledrydep.PM25) {
	return simpledrydep.SOx{igS}, simpledrydep.NH3{igNH}, simpledrydep.NOx{igNO}, simpledrydep.VOC{igOrg}, simpledrydep.PM25{ipOrg, iPM2_5, ipNH, ipS, ipNO}
}

func TestDryDeposition(t *testing.T) {
	cfg, ctmdata, pop, popIndices, mr, mortIndices := inmap.VarGridTestData()
	emis := inmap.NewEmissions()

	mutator, err := inmap.PopulationMutator(cfg, popIndices)
	if err != nil {
		t.Error(err)
	}
	var m simplechem.Mechanism
	d := &inmap.InMAP{
		InitFuncs: []inmap.DomainManipulator{
			cfg.RegularGrid(ctmdata, pop, popIndices, mr, mortIndices, emis, m),
			cfg.MutateGrid(mutator, ctmdata, pop, mr, emis, m, nil),
			inmap.SetTimestepCFL(),
		},
		RunFuncs: []inmap.DomainManipulator{
			inmap.Calculations(inmap.AddEmissionsFlux()),
			inmap.Calculations(simpledrydep.DryDeposition(simpleDryDepIndices)),
			inmap.SteadyStateConvergenceCheck(1, cfg.PopGridColumn, m, nil),
		},
	}
	if err := d.Init(); err != nil {
		t.Error(err)
	}
	for _, c := range d.Cells() {
		for i := range c.Ci {
			c.Cf[i] = 1 // set concentrations to 1
		}
	}
	if err := d.Run(); err != nil {
		t.Error(err)
	}

	for _, c := range d.Cells() {
		for ii, cc := range c.Cf {
			if c.Layer == 0 {
				if cc >= 1 || cc <= 0.98 {
					t.Errorf("ground-level cell %v pollutant %d should equal be between 0.98 and 1 but is %g", c, ii, cc)
				}
			} else if cc != 1 {
				t.Errorf("above-ground cell %v pollutant %d should equal 1 but equals %g", c, ii, cc)
			}
		}
	}
}
