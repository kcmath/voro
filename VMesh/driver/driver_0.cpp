#include <iostream>
#include <string>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <memory>
#include <vtkXMLDataSetWriter.h>
#include <vtkIntArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "voro++.hh"

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

// functions to create cell data objects using vtkSmartPointer
vtkSmartPointer<vtkIntArray> createCellAttributeInt(int nComponents, int nTuples, const char* attName) {
   vtkSmartPointer<vtkIntArray> cellAttribute = vtkSmartPointer<vtkIntArray>::New();
   cellAttribute->SetNumberOfComponents(nComponents);
   cellAttribute->SetNumberOfTuples(nTuples);
   cellAttribute->SetName(attName);
   return cellAttribute;
}

vtkSmartPointer<vtkDoubleArray> createCellAttributeDouble(int nComponents, int nTuples, const char* attName) {
   vtkSmartPointer<vtkDoubleArray> cellAttribute = vtkSmartPointer<vtkDoubleArray>::New();
   cellAttribute->SetNumberOfComponents(nComponents);
   cellAttribute->SetNumberOfTuples(nTuples);
   cellAttribute->SetName(attName);
   return cellAttribute;
}


bool is_big_endian(void) {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}


int main(int argc, const char * argv[])
{

// Set up constants for the container geometry
const double x_min=-1,x_max=1;
const double y_min=-1,y_max=1;
const double z_min=-1,z_max=1;
const double cvol=(x_max-x_min)*(y_max-y_min)*(z_max-z_min);

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=20;

	int i;
	double x,y,z;


	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Randomly add particles into the container
	for(i=0;i<particles;i++) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		con.put(i,x,y,z);
	}


      // create unstructured grid and points (i.e. vertices)
      vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

      // create cell attributes
      vtkSmartPointer<vtkIntArray> cellID = createCellAttributeInt(1, particles, "cellID");
      vtkSmartPointer<vtkIntArray> cellFaces = createCellAttributeInt(1, particles, "cellFaces");
      vtkSmartPointer<vtkDoubleArray> cellVolume = createCellAttributeDouble(1, particles, "cellVolume");
      vtkSmartPointer<vtkDoubleArray> cellSurfaceArea = createCellAttributeDouble(1, particles, "cellSurfaceArea");
      
      // total number of vertices in the container with duplicates
      int containerNumberOfVerticesWithDups = 0;

      // initialize variables used in cell loop
      voro::voronoicell_neighbor c;        // create a Voronoi cell with neighbor information
      std::vector<int> neigh;        // neighbors' information
      std::vector<int> f_vert;       // vertices indexes
      std::vector<int> f_order;      // number of vertices per face
      std::vector<double> v;

      // create container loop object
      voro::c_loop_all cellLoop(con);      
      // loop through cells
      if( cellLoop.start() ) {
         int counter = 0;
         do {
            if( con.compute_cell(c,cellLoop) ) {  // get Voronoi cell information

               // gather information about the computed Voronoi cell
               cellLoop.pos(x,y,z);
               c.neighbors(neigh);
               c.face_vertices(f_vert);
               c.vertices(x,y,z,v);

               // loop vertices and store their position
               for( unsigned int i=0 ; i<v.size() ; i+=3 ) {
                  points->InsertNextPoint(v[i], v[i+1], v[i+2]);
               }

               // vtk faces
               // [numberOfCellFaces, (numberOfPointsOfFace0, pointId0, pointId1, … ),
               //                     (numberOfPointsOfFace1, pointId0, pointId1, …), … ].
               // create faces ID list
               vtkSmartPointer<vtkIdList> vtkFaces = vtkSmartPointer<vtkIdList>::New();
               vtkFaces->InsertNextId(neigh.size());

               // number of vertices in current cell
               int numberOfVertices = v.size()/3;

               // update total number of vertices in container
               containerNumberOfVerticesWithDups += numberOfVertices;

               // update starting index for the current cell in the container
               int containerVertexStartIndex = containerNumberOfVerticesWithDups - numberOfVertices;

               // loop over all faces of the Voronoi cell and populate vtkFaces with
               // numberOfVerticesPerFace and their vertex indices
               int j,k=0;
               int numberOfVerticesPerFace;
               while( (unsigned int)k<f_vert.size() ) {
                  numberOfVerticesPerFace = f_vert[k++];
                  vtkFaces->InsertNextId(numberOfVerticesPerFace);  // number of vertices in 1 face

                  j = k+numberOfVerticesPerFace;
                  while( k<j ) {
                     int containerIndex = f_vert[k++] + containerVertexStartIndex;
                     vtkFaces->InsertNextId(containerIndex);
                  } // end single face loop
               } // end vertices loop

               // add cell to unstructure grid
               uGrid->InsertNextCell(VTK_POLYHEDRON,vtkFaces);

               // add attributes to cell
               cellID->InsertValue(counter, counter);
               cellVolume->InsertValue(counter, c.volume());
               cellFaces->InsertValue(counter, neigh.size());
               cellSurfaceArea->InsertValue(counter, c.surface_area());

               counter++;
            } // end neighbor compute if
         } while(cellLoop.inc()); // end do loop

      } // end cell loop if

      // add attributes to cellData
      uGrid->GetCellData()->AddArray(cellID);
      uGrid->GetCellData()->AddArray(cellVolume);
      uGrid->GetCellData()->AddArray(cellFaces);
      uGrid->GetCellData()->AddArray(cellSurfaceArea);

      // add points to unstructured grid
      uGrid->SetPoints(points);

      // create .vtu file name
      char filename[256];
      sprintf(filename,"random.vtu");

      // output unstructured grid
      std::cout << "Writing .vtu file..." << std::endl;
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtkWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      vtkWriter->SetInputData(uGrid);
      vtkWriter->SetFileName(filename);
      //vtkWriter->SetDataModeToAscii();
      vtkWriter->SetDataModeToBinary();  // much smaller files and faster
      vtkWriter->Update();
      

	// Sum up the volumes, and check that this matches the container volume
	double vvol=con.sum_cell_volumes();
	printf("Container volume : %g\n"
	       "Voronoi volume   : %g\n"
	       "Difference       : %g\n",cvol,vvol,vvol-cvol);

	// Output the particle positions in gnuplot format
	con.draw_particles("random_points_p.gnu");

	// Output the Voronoi cells in gnuplot format
	con.draw_cells_gnuplot("random_points_v.gnu");

    return 0;
}

