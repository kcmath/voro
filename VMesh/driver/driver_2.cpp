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
#include <Eigen/Dense>

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

	voro::voronoicell v;
   int n_faces;                             // number of faces
   double x,y,z;                            // centroid position
   std::vector<double> vertices_positions;  // vertices' positions
   std::vector<int> f_vert;                 // vertices indexes

	// Initialize the Voronoi cell to be a cube of side length 2, centered
	// on the origin
	v.init(-1,1,-1,1,-1,1);

	// Remove one edge of the cell with a single plane cut
	v.plane(1,1,0,2);

	// Output the Voronoi cell to a file in gnuplot format
	v.draw_gnuplot(0,0,0,"simple_cell.gnu");

	// Output vertex-based statistics
	printf("Total vertices      : %d\n",v.p);
	printf("Vertex positions    : ");v.output_vertices();puts("");
	printf("Vertex orders       : ");v.output_vertex_orders();puts("");
	printf("Max rad. sq. vertex : %g\n\n",0.25*v.max_radius_squared());

	// Output edge-based statistics
	printf("Total edges         : %d\n",v.number_of_edges());
	printf("Total edge distance : %g\n",v.total_edge_distance());
	printf("Face perimeters     : ");v.output_face_perimeters();puts("\n");

	// Output face-based statistics
	printf("Total faces         : %d\n",v.number_of_faces());
	printf("Surface area        : %g\n",v.surface_area());
	printf("Face freq. table    : ");v.output_face_freq_table();puts("");
	printf("Face orders         : ");v.output_face_orders();puts("");
	printf("Face areas          : ");v.output_face_areas();puts("");
	printf("Face normals        : ");v.output_normals();puts("");
	printf("Face vertices       : ");v.output_face_vertices();puts("\n");

	// Output volume-based statistics
	v.centroid(x,y,z);
	printf("Volume              : %g\n"
	       "Centroid vector     : (%g,%g,%g)\n",v.volume(),x,y,z);

   // create unstructured grid and points (i.e. vertices)
   vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

   // gather information about the computed Voronoi cell
   n_faces = v.number_of_faces();
   v.face_vertices(f_vert);
   v.vertices(vertices_positions);

   // loop vertices and store their position
   for( unsigned int i=0 ; i<vertices_positions.size() ; i+=3 ) {
      points->InsertNextPoint(vertices_positions[i], vertices_positions[i+1], vertices_positions[i+2]);
      }
      
   // vtk faces
   // [numberOfCellFaces, (numberOfPointsOfFace0, pointId0, pointId1, … ),
   //                     (numberOfPointsOfFace1, pointId0, pointId1, …), … ].
   // create faces ID list
   vtkSmartPointer<vtkIdList> vtkFaces = vtkSmartPointer<vtkIdList>::New();
   vtkFaces->InsertNextId(n_faces);

   // loop over all faces of the Voronoi cell and populate vtkFaces with
   // numberOfVerticesPerFace and their vertex indices
   int j,k=0;
   int numberOfVerticesPerFace;
   while( (unsigned int)k<f_vert.size() ) {
      numberOfVerticesPerFace = f_vert[k++];
      vtkFaces->InsertNextId(numberOfVerticesPerFace);  // number of vertices in 1 face

      j = k+numberOfVerticesPerFace;
      while( k<j ) {
         int containerIndex = f_vert[k++];
         vtkFaces->InsertNextId(containerIndex);
      } // end single face loop
   } // end vertices loop

   // add cell to unstructure grid
   uGrid->InsertNextCell(VTK_POLYHEDRON,vtkFaces);

   // add points to unstructured grid
   uGrid->SetPoints(points);

   // output unstructured grid
   std::cout << "Writing .vtu file..." << std::endl;
   vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtkWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
   vtkWriter->SetInputData(uGrid);
   vtkWriter->SetFileName("cell_statistics_vtk.vtu");
   vtkWriter->SetDataModeToAscii();
   //vtkWriter->SetDataModeToBinary();
   vtkWriter->Update();
         
    return 0;
}

