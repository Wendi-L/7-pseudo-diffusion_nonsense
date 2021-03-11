/*****************************************************************************	
* Multiscale Universal Interface Code Coupling Library Demo 7                *	
*                                                                            *	
* Copyright (C) 2019 W. Liu                                                  *	
*                                                                            *	
* This software is jointly licensed under the Apache License, Version 2.0    *	
* and the GNU General Public License version 3, you may use it according     *	
* to either.                                                                 *	
*                                                                            *	
* ** Apache License, version 2.0 **                                          *	
*                                                                            *	
* Licensed under the Apache License, Version 2.0 (the "License");            *	
* you may not use this file except in compliance with the License.           *	
* You may obtain a copy of the License at                                    *	
*                                                                            *	
* http://www.apache.org/licenses/LICENSE-2.0                                 *	
*                                                                            *	
* Unless required by applicable law or agreed to in writing, software        *	
* distributed under the License is distributed on an "AS IS" BASIS,          *	
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *	
* See the License for the specific language governing permissions and        *	
* limitations under the License.                                             *	
*                                                                            *	
* ** GNU General Public License, version 3 **                                *	
*                                                                            *	
* This program is free software: you can redistribute it and/or modify       *	
* it under the terms of the GNU General Public License as published by       *	
* the Free Software Foundation, either version 3 of the License, or          *	
* (at your option) any later version.                                        *	
*                                                                            *	
* This program is distributed in the hope that it will be useful,            *	
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *	
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *	
* GNU General Public License for more details.                               *	
*                                                                            *	
* You should have received a copy of the GNU General Public License          *	
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *	
******************************************************************************/	

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

/// Include MUI header file and configure file 
#include "mui.h"
#include "demo7_config.h"

int main(int argc, char ** argv) {

    /// Create results folder
    mkdir("coupling_results", 0777);

    /// Create rbf matrix folder
    mkdir("rbfFineMatrix", 0777);

    /// Declare MPI common world with the scope of MUI
	MPI_Comm  world = mui::mpi_split_by_app();

    /// Declare MPI ranks and rank size
    int rank, size;
    MPI_Comm_rank( world, &rank );
    MPI_Comm_size( world, &size );

    /// Define the name of MUI domain
    std::string domain = "fineDomain";

    /// Define the name of MUI interfaces 
    std::vector<std::string> interfaces;
    interfaces.emplace_back( "interface2D01" );
    interfaces.emplace_back( "interface2D02" );

    /// Declare MUI objects using MUI configure file
    auto ifs = mui::create_uniface<mui::demo7_config>( domain, interfaces );

    /// Define the name of push/fetch values
	const char* name_push   = "fineField";			
	const char* name_fetch  = "coarseField";

    /// Define the forget steps of MUI to reduce the memory
    int         forgetSteps = 5;

    /// Define parameters of the RBF sampler
    /// Define the search radius of the RBF sampler
    /// The search radius should not set to a very large value so that to ensure a good convergence
    double      rSampler    = 	1000;
	bool 		conservative=	false; 
	double 		cutoff		=	1e-9;
	bool 		polynomial	=	true;
	bool 		smoothFunc	=	false;
	bool 		readMatrix	=	false;
    std::string fileAddress("rbfFineMatrix");

	/// Setup time steps
    constexpr static int    steps           = 200;

    /// Setup output interval
    constexpr static int    outputInterval  = 200;
	
    constexpr static int    Nt2      = 100;  /// total number of points
    double points2[Nt2][3], points3[Nt2][3];
	
    /// Store point coordinates
 	for ( int k = 0; k < Nt2; ++k ) {
                points2[k][0] = k;
                points2[k][1] = k;
                points2[k][2] = k;

                points3[k][0] = k;
                points3[k][1] = k;
                points3[k][2] = k;
    }

  	double scalar_field2[Nt2], scalar_field3[Nt2];

 	for ( int k = 0; k < Nt2; ++k ) {
                scalar_field2[k] = 100.0;
    }

    /// Declare std::vector to store mui::point2d
    std::vector<mui::point3d> point3dvec;

    /// Store mui::point2d that located in the fetch interface
 	for ( int k = 0; k < Nt2; ++k ) {
                    mui::point3d ptf( points3[k][0], points3[k][1], points3[k][2]);
                    point3dvec.push_back(ptf);
    }

    /// Define and announce MUI send/receive span
    mui::geometry::box<mui::demo7_config> send_region( {0, 0, 0}, {99, 99, 99} );
    mui::geometry::box<mui::demo7_config> recv_region( {0, 0, 0}, {99, 99, 99} );
    ifs[0]->announce_send_span( 0, steps, send_region );
    ifs[1]->announce_recv_span( 0, steps, recv_region );

	/// Define spatial and temporal samplers
    mui::sampler_rbf<mui::demo7_config> spatial_sampler(rSampler,point3dvec,conservative,cutoff,polynomial,fileAddress,readMatrix);
	//mui::sampler_pseudo_nearest_neighbor<mui::demo7_config> spatial_sampler(rSampler);
    mui::chrono_sampler_exact<mui::demo7_config> chrono_sampler;

	/// Commit ZERO step of MUI
	ifs[0]->commit(0);

    /// Output the initial pseudo scalar field
    std::ofstream outputFileLeft;
    std::string filenameL = "coupling_results/scalar_field_left_fine_0.csv";
    outputFileLeft.open (filenameL);
    outputFileLeft << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";
 	for ( int k = 0; k < Nt2; ++k ) {
                outputFileLeft << points2[k][0] << "," <<points2[k][1]<< "," << points2[k][2]<< "," << scalar_field2[k] << ", \n";
            }

    outputFileLeft.close();

    std::ofstream outputFileRight;
    std::string filenameR = "coupling_results/scalar_field_right_fine_0.csv";
    outputFileRight.open (filenameR);
    outputFileRight << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n"; 
 	for ( int k = 0; k < Nt2; ++k ) {
                outputFileRight << points3[k][0] << "," <<points3[k][1]<< "," << points3[k][2]<< "," << scalar_field3[k] << ", \n";

    }
    outputFileRight.close();

    /// Define output files for boundary integrations
    std::ofstream outputIntegration;
    std::string outputIntegrationName = "coupling_results/faceIntegrationD1_D3.txt";
    outputIntegration.open (outputIntegrationName);
    outputIntegration << "\"t\",\"intFaceLD1\",\"intFaceRD1\",\"intFaceLD3\",\"intFaceRD3\"\n";
    outputIntegration.close();

	/// Begin time loops
    for ( int t = 1; t <= steps; ++t ) {

		printf("\n");
		printf("{Fine Domain} %d Step \n", t );

        /// Loop over points of Domain 1
        for ( int k = 0; k < Nt2; ++k ) {
			mui::point3d locp( points2[k][0], points2[k][1], points2[k][2] );
			/// push data to the other solver
			ifs[0]->push( name_push, locp, scalar_field2[k] );
        }

        /// Commit 't' step of MUI
        int sent = ifs[0]->commit( t );

        /// Forget data from Zero step to 't - forgetSteps' step of MUI to save memory
        if ( (t - forgetSteps) > 0 ) {
            ifs[0]->forget( t - forgetSteps );
        }

        /// Loop over points of Domain 3
        for ( int k = 0; k < Nt2; ++k ) {
			mui::point3d locf( points3[k][0], points3[k][1], points3[k][2]);
			/// Fetch data from the other solver
			scalar_field3[k] = ifs[1]->fetch( name_fetch,
													locf,
													t,
													spatial_sampler,
													chrono_sampler);
		}

        /// Output the pseudo scalar field and the boundary integrations
        if ((t % outputInterval) == 0) {

            std::ofstream outputFileLeft;
            std::string filenameL = "coupling_results/scalar_field_left_fine_" + std::to_string(t) + ".csv";
            outputFileLeft.open (filenameL);
            outputFileLeft << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";

            for ( int k = 0; k < Nt2; ++k ) {
                        outputFileLeft << points2[k][0] << "," <<points2[k][1]<< "," << points2[k][2]<< "," << scalar_field2[k] << ", \n";
            }
            outputFileLeft.close();

            std::ofstream outputFileRight;
            std::string filenameR = "coupling_results/scalar_field_right_fine_" + std::to_string(t) + ".csv";
            outputFileRight.open (filenameR);
            outputFileRight << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";
            for ( int k = 0; k < Nt2; ++k ) {
                        outputFileRight << points3[k][0] << "," <<points3[k][1]<< "," << points3[k][2]<< "," << scalar_field3[k] << ", \n";
            }
            outputFileRight.close();

        }
	}

    return 0;
}