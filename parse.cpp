#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vtkVertex.h>

struct Point {
    int id;
    float x, y, z;

    Point(int id, float x, float y, float z) : id(id), x(x), y(y), z(z) {}
};

int main() {
    std::ifstream inFile("yourfile.dat");
    if (!inFile) {
        std::cerr << "Unable to open file";
        exit(1);
    }

    std::vector<Point> points;
    std::string line;
    int id = 0; // Initialize an ID counter

    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        char brace, comma;
        float x, y, z;

        while (ss >> brace >> x >> comma >> y >> comma >> z >> brace) {
            points.emplace_back(id++, x, y, z); // Add new Point to vector with a unique ID
        }
    }

    inFile.close();

    // Print out the points to verify
    for (const auto& point : points) {
        std::cout << "ID: " << point.id << " - Point: {" << point.x << ", " << point.y << ", " << point.z << "}" << std::endl;
    }

vtkSmartPointer<vtkPoints> ReadData(const std::string &fileName) {
    std::ifstream inFile(fileName);
    if (!inFile) {
        std::cerr << "Unable to open file";
        exit(1);
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    std::string line;

    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        char brace, comma;
        double x, y, z;

        while (ss >> brace >> x >> comma >> y >> comma >> z >> brace) {
            points->InsertNextPoint(x, y, z);
        }
    }

    inFile.close();
    return points;
}

int main() {
    // Replace with the path to your .dat file
    std::string inputFileName = "yourfile.dat";

    // Read the points from the file
    vtkSmartPointer<vtkPoints> points = ReadData(inputFileName);

    // Create a polydata object and set the points
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);

    // Write the polydata to a .vtp file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output.vtp");
    writer->SetInputData(polydata);

    // Optional: Set the data mode (Ascii, Binary, etc.)
    writer->SetDataModeToBinary();

    // Write the file
    writer->Write();

    return 0;
}

vtkSmartPointer<vtkPoints> ReadData(const std::string &fileName) {
    std::ifstream inFile(fileName);
    if (!inFile) {
        std::cerr << "Unable to open file";
        exit(1);
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    std::string line;

    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        char brace, comma;
        double x, y, z;

        while (ss >> brace >> x >> comma >> y >> comma >> z >> brace) {
            points->InsertNextPoint(x, y, z);
        }
    }

    inFile.close();
    return points;
}

int main() {
    // Replace with the path to your .dat file
    std::string inputFileName = "yourfile.dat";

    // Read the points from the file
    vtkSmartPointer<vtkPoints> points = ReadData(inputFileName);

    // Create an unstructured grid object and set the points
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);

    // Add vertices to the grid for each point
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
        vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
        vertex->GetPointIds()->SetId(0, i);
        unstructuredGrid->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
    }

    // Write the unstructured grid to a .vtu file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("output.vtu");
    writer->SetInputData(unstructuredGrid);

    // Optional: Set the data mode (Ascii, Binary, etc.)
    writer->SetDataModeToBinary();

    // Write the file
    writer->Write();

    return 0;
}

    
