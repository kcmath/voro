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
#include <VoroMesh.hpp>
#include <Eigen/Dense>

// Major and minor torus radii
const double arad=9,brad=3.5;

// The outer radius of the torus, that determines how big the container should
// be
const double crad=arad+brad;

// Set up constants for the container geometry
const double x_min=-crad-0.5,x_max=crad+0.5;
const double y_min=-crad-0.5,y_max=crad+0.5;
const double z_min=-brad-0.5,z_max=brad+0.5;

// Set the computational grid size
const int n_x=10,n_y=10,n_z=3;

// This class creates a custom toroidal wall object that is centered on the
// origin and is aligned with the xy plane. It is derived from the pure virtual
// "wall" class. The "wall" class contains virtual functions for cutting the
// Voronoi cell in response to a wall, and for telling whether a given point is
// inside the wall or not. In this derived class, specific implementations of
// these functions are given.
class wall_torus : public voro::wall {
	public:

		// The wall constructor initializes constants for the major and
		// minor axes of the torus. It also initializes the wall ID
		// number that is used when the plane cuts are made. This is
		// only tracked with the voronoicell_neighbor class and is
		// ignored otherwise. It can be omitted, and then an arbitrary
		// value of -99 is used.
		wall_torus(double imjr,double imnr,int iw_id=-99)
			: w_id(iw_id), mjr(imjr), mnr(imnr) {};

		// This returns true if a given vector is inside the torus, and
		// false if it is outside. For the current example, this
		// routine is not needed, but in general it would be, for use
		// with the point_inside() routine in the container class.
		bool point_inside(double x,double y,double z) {
			double temp=sqrt(x*x+y*y)-mjr;
			return temp*temp+z*z<mnr*mnr;
		}

		// This template takes a reference to a voronoicell or
		// voronoicell_neighbor object for a particle at a vector
		// (x,y,z), and makes a plane cut to to the object to account
		// for the toroidal wall
		template<class vc_class>
		inline bool cut_cell_base(vc_class &c,double x,double y,double z) {
			double orad=sqrt(x*x+y*y);
			double odis=orad-mjr;
			double ot=odis*odis+z*z;

			// Unless the particle is within 1% of the minor
			// radius, then a plane cut is made
			if(ot>0.01*mnr*mnr) {
				ot=2*mnr/sqrt(ot)-2;
				z*=ot;
				odis*=ot/orad;
				x*=odis;
				y*=odis;
				return c.nplane(x,y,z,w_id);
			}
			return true;
		}

		// These virtual functions are called during the cell
		// computation in the container class. They call instances of
		// the template given above.
		bool cut_cell(voro::voronoicell &c,double x,
				double y,double z) {return cut_cell_base(c,x,y,z);}
		bool cut_cell(voro::voronoicell_neighbor &c,double x,
				double y,double z) {return cut_cell_base(c,x,y,z);}
	private:
		// The ID number associated with the wall
		const int w_id;
		// The major radius of the torus
		const double mjr;
		// The minor radius of the torus
		const double mnr;
};

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
double x,y,z;   // particle position
// Set the number of particles that are going to be randomly introduced
const int particles=2393;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Add the custom toroidal wall to the container
	wall_torus tor(arad,brad);
	con.add_wall(tor);

	// Import the particles from a file
	con.import("pack_torus");

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
      sprintf(filename,"torus.vtu");

      // output unstructured grid
      std::cout << "Writing .vtu file..." << std::endl;
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtkWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      vtkWriter->SetInputData(uGrid);
      vtkWriter->SetFileName(filename);
      vtkWriter->SetDataModeToAscii();
      //vtkWriter->SetDataModeToBinary();  // much smaller files and faster
      vtkWriter->Update();      
	
    return 0;
}

