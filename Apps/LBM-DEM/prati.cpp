/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice
 *    this list of conditions and the following disclaimer in the documentation
 *    and or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
*/

/*! @brief   Implement member functions for PRATI scheme
 * @author C. Tsigginos
 * @details Implementing functions for PRATI.
 */

#include "prati.h"
#include "flowfield.h"
#include "flowfield_host_device.h"
#include "scheme.h"
#include <vector>
#include <map>
#include "ops_seq_v2.h"
#include "model_kernel.inc"
#include "prati.inc"
Prati::Prati(Component componentUser, int spacedim, bool owned = false,
		int porosModel = 0, Real gammaUser = 0.0, int nelem) :
		FsiBase(componentUser, spacedim, owned, porosModel),
		Fd{"FdPrati"} {

	noElem = nelem;

	if ((porosModel!= Mode_Spherical) || (porosModel != Mode_Grid)) {
		ops_printf("ERROR: Implemented porosity model not consistent with PRATI scheme\n");
		exit(EXIT_FAILURE);
	}

	if (porosModel == Mode_Spherical) {
		poros = new PorosSpherical();
	}
	else if (porosModel == Mode_Grid)
		poros = new PorosGrid();


}

void Prati::DefineVariables(SizeType timestep) {


	int size = noElem * spaceDim;

	Fd.SetDataDim(size);
	if (timestep == 0) {
		Fd. CreateFieldFromScratch(g_Block());
		poros.DefineVariables(noElem);
	}

	else {
		Fd.CreateFieldFromFile(CaseName(), g_Block(), timestep);
		poros.DefineVariables(noElem, timestep);
	}

}

void Prati::ModelCollision() {

	int size = noElem * spaceDim;



	for (const auto& idBlock : g_Block()) {
		const Block& block{idBlock.second};
		std::vector<int> iterRng;
		iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());

		const Block& block{idBlock.second};
		Real tauRef = compo.tauRef();
		const Real* pdt {pTimestep()};
		const int blockIndex{block.ID()};
		ops_par_loop(KerCollisionPrati3D,"KerCollisionPrati3D", block.Get(),
					 spaceDim, iterRng.data(),
					 ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
							 	 "double", OPS_WRITE),
					 ops_arg_dat(Fd[blockIndex], size, LOCALSTENCIL,
							 	 "double", OPS_WRITE),
					 ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
							 	 "double", OPS_READ),
					 ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
                                 LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
							 	 LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_Rho).id)
								 .at(blockIndex), 1, LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
								 1, LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
								 1, LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIndex),
								 1, LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(poros->id[blockIndex], noElem, LOCALSTENCIL,
							 "int", OPS_READ),
					 ops_arg_dat(poros->sfParticle[blockIndex], noElem, LOCALSTENCIL,
							 	 "double", OPS_READ),
					 ops_arg_dat(poros->vParts[blockIndex], size, LOCALSTENCIL),
					 	 	 	 "double", OPS_READ),
					 ops_arg_gbl(pdt, 1, "double", OPS_READ),
					 ops_arg_gbl(&tauRef, 1, "int", OPS_READ),
					 ops_arg_gbl(&forceFlag, 1, "int", OPS_READ),
					 ops_arg_gbl(force, spaceDim, "double", OPS_READ),
					 ops_arg_gbl(&gamma, 1, "double", OPS_READ),
					 ops_arg_gbl(&noElem, 1, "int", OPS_READ),
					 ops_arg_gbl(&spaceDim, 1, "int", OPS_READ)
					 ops_arg_gbl(compo.index, 2, "int", OPS_READ));
	}



}

void Prati::InitializeVariables() {

	int size = noElem * spaceDim;

	poros->InitializeVariables();

	for (const auto& idBlock : g_Block()) {
		const Block& block{idBlock.second};
		std::vector<int> iterRng;
		iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());

		const int blockIndex{block.ID()};
		ops_par_loop(KerInitialize,"KerInitialize", block.Get(),
					 spaceDim, iterRng.data(),
					 ops_arg_dat(Fd[blockIndex], noElem, LOCALSTENCIL,
							 	 "double", OPS_WRITE),
					 ops_arg_gbl(&size, 1, "int", OPS_READ));

	}
}

void Prati::PostVelocityCalculation() {

	Real tauRef = compo.tauRef();
	const Real* pdt {pTimestep()};

	int size = noElem * spaceDim;
	for (const auto& idBlock : g_Block()) {
		const Block& block(idBlock.second());
		 std::vector<int> iterRng;
		 iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
		 const int blockIndex{block.ID()};
		 for (auto& macroVar : compo.macroVars) {
			 const int varId{macroVar.second.id};
			 const VariableTypes varType{macroVar.first};
			 switch (varType) {
			 	 case Variable_U:
			 		 ops_par_loop(KerPratiUpdateU3D, "KerPratiUpdateU3D", block.Get(),
			 				 	  spaceDim, iterRng.data,
								  ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
	                                          1, LOCALSTENCIL, "double", OPS_RW),
								  ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
											  LOCALSTENCIL, "int", OPS_READ),
								  ops_arg_dat(poros->id[blockIndex], noElem,
										  	  LOCALSTENCIL, "int", OPS_READ),
								  ops_arg_dat(poros->poros[blockIndex], noElem,
										  	  LOCALSTENCIL, "double", OPS_READ),
								  ops_arg_dat(poros->vParts[blockIndex], size,
										  	  LOCALSTENCIL, "double", OPS_READ),
								  ops_arg_gbl(pdt, 1, "double", OPS_READ),
								  ops_arg_gbl(&tauRef, 1, "double", OPS_READ),
								  ops_arg_gbl(&gamma, 1, "double", OPS_READ),
								  ops_arg_gbl(&forceFlag, 1, "int", OPS_READ),
								  ops_arg_gbl(force, spaceDim, "double", OPS_READ),
								  ops_arg_gbl(&noElem, 1, "int", OPS_READ),
								  ops_arg_gbl(&spaceDim, 1, "int", OPS_READ));
			 		 break;
			 	 case Variable_V:
			 		 ops_par_loop(KerPratiUpdateV3D, "KerPratiUpdateV3D", block.Get(),
			 					  spaceDim, iterRng.data,
			 					  ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
			 			                      1, LOCALSTENCIL, "double", OPS_RW),
								  ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
											  LOCALSTENCIL, "int", OPS_READ),
			 					  ops_arg_dat(poros->id[blockIndex], noElem,
			 								  LOCALSTENCIL, "int", OPS_READ),
			 				      ops_arg_dat(poros->poros[blockIndex], noElem,
			 						          LOCALSTENCIL, "double", OPS_READ),
			 				      ops_arg_dat(poros->vParts[blockIndex], size,
			 								  LOCALSTENCIL, "double", OPS_READ),
								  ops_arg_gbl(pdt, 1, "double", OPS_READ),
								  ops_arg_gbl(&tauRef, 1, "double", OPS_READ),
								  ops_arg_gbl(&gamma, 1, "double", OPS_READ),
			 					  ops_arg_gbl(&forceFlag, 1, "int", OPS_READ),
			 					  ops_arg_gbl(force, spaceDim, "double", OPS_READ),
			 					  ops_arg_gbl(&noElem, 1, "int", OPS_READ),
			 					  ops_arg_gbl(&spaceDim, 1, "int", OPS_READ));
			 		 break;
			 	 case Variable_W:
			 		 ops_par_loop(KerPratiUpdateW3D, "KerPratiUpdateW3D", block.Get(),
			 					  spaceDim, iterRng.data,
			 					  ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
			 			                      1, LOCALSTENCIL, "double", OPS_RW),
								  ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
											  LOCALSTENCIL, "int", OPS_READ),
			 					  ops_arg_dat(poros->id[blockIndex], noElem,
			 								  LOCALSTENCIL, "int", OPS_READ),
			 				      ops_arg_dat(poros->poros[blockIndex], noElem,
			 						          LOCALSTENCIL, "double", OPS_READ),
			 				      ops_arg_dat(poros->vParts[blockIndex], size,
			 								  LOCALSTENCIL, "double", OPS_READ),
								  ops_arg_gbl(pdt, 1, "double", OPS_READ),
								  ops_arg_gbl(&tauRef, 1, "double", OPS_READ),
								  ops_arg_gbl(&gamma, 1, "double", OPS_READ),
			 					  ops_arg_gbl(&forceFlag, 1, "int", OPS_READ),
			 					  ops_arg_gbl(force, spaceDim, "double", OPS_READ),
			 					  ops_arg_gbl(&noElem, 1, "int", OPS_READ),
			 					  ops_arg_gbl(&spaceDim, 1, "int", OPS_READ));
			 		 break;
			 	 default:
			 		 break;

			 }
		 }
	}

}
