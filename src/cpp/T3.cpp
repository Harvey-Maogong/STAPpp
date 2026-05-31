/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     T3 Element Extension                                                  */
/*****************************************************************************/

#include "T3.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CT3::CT3()
{
    NEN_ = 3;			// Each element has 3 nodes
    nodes_ = new CNode * [NEN_];

    ND_ = 6;			// 3 nodes * 2 DOFs (ux, uy) for plane problems
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//	Desconstructor
CT3::~CT3()
{
}

//	Read element data from stream Input
bool CT3::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;		// Material property set number
    unsigned int N1, N2, N3;	// Node numbers

    Input >> N1 >> N2 >> N3 >> MSet;
    ElementMaterial_ = dynamic_cast <CT3Material*> (MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];

    return true;
}

//	Write element data to stream
void CT3::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
        << setw(9) << nodes_[1]->NodeNumber
        << setw(9) << nodes_[2]->NodeNumber
        << setw(12) << ElementMaterial_->nset << endl;
}

//	Generate location matrix: map element DOFs to global equation numbers
void CT3::GenerateLocationMatrix()
{
    // Plane T3 element uses only x and y DOFs of each node
    LocationMatrix_[0] = nodes_[0]->bcode[0];	// Node 1, x
    LocationMatrix_[1] = nodes_[0]->bcode[1];	// Node 1, y
    LocationMatrix_[2] = nodes_[1]->bcode[0];	// Node 2, x
    LocationMatrix_[3] = nodes_[1]->bcode[1];	// Node 2, y
    LocationMatrix_[4] = nodes_[2]->bcode[0];	// Node 3, x
    LocationMatrix_[5] = nodes_[2]->bcode[1];	// Node 3, y
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by column starting from the diagonal element
void CT3::ElementStiffness(double* Matrix)
{
    clear(Matrix, 21);

    CT3Material* material = dynamic_cast <CT3Material*> (ElementMaterial_);

    //	Node coordinates
    double x1 = nodes_[0]->XYZ[0];
    double y1 = nodes_[0]->XYZ[1];
    double x2 = nodes_[1]->XYZ[0];
    double y2 = nodes_[1]->XYZ[1];
    double x3 = nodes_[2]->XYZ[0];
    double y3 = nodes_[2]->XYZ[1];

    //	Shape function derivatives (area coordinates coefficients)
    double b1 = y2 - y3;
    double b2 = y3 - y1;
    double b3 = y1 - y2;
    double c1 = x3 - x2;
    double c2 = x1 - x3;
    double c3 = x2 - x1;

    double Area2 = b1 * c2 - b2 * c1;	// 2 * Area

    if (Area2 <= 0.0)
    {
        cerr << "*** Error *** T3 element has non-positive area." << endl;
        return;
    }

    double Area = Area2 * 0.5;
    double inv2A = 1.0 / Area2;

    //	B matrix (3 x 6) - constant strain
    double B[3][6];
    B[0][0] = b1 * inv2A;  B[0][1] = 0.0;        B[0][2] = b2 * inv2A;  B[0][3] = 0.0;        B[0][4] = b3 * inv2A;  B[0][5] = 0.0;
    B[1][0] = 0.0;         B[1][1] = c1 * inv2A;  B[1][2] = 0.0;         B[1][3] = c2 * inv2A;  B[1][4] = 0.0;         B[1][5] = c3 * inv2A;
    B[2][0] = c1 * inv2A;  B[2][1] = b1 * inv2A;  B[2][2] = c2 * inv2A;  B[2][3] = b2 * inv2A;  B[2][4] = c3 * inv2A;  B[2][5] = b3 * inv2A;

    //	D matrix (3 x 3) - plane stress / plane strain
    double D[3][3];
    double E = material->E;
    double nu = material->nu;
    double factor;

    if (material->type == 1)	// Plane stress
    {
        factor = E / (1.0 - nu * nu);
        D[0][0] = factor;        D[0][1] = factor * nu;  D[0][2] = 0.0;
        D[1][0] = factor * nu;   D[1][1] = factor;       D[1][2] = 0.0;
        D[2][0] = 0.0;           D[2][1] = 0.0;          D[2][2] = factor * 0.5 * (1.0 - nu);
    }
    else	// Plane strain
    {
        factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        D[0][0] = factor * (1.0 - nu);  D[0][1] = factor * nu;          D[0][2] = 0.0;
        D[1][0] = factor * nu;          D[1][1] = factor * (1.0 - nu);  D[1][2] = 0.0;
        D[2][0] = 0.0;                  D[2][1] = 0.0;                  D[2][2] = factor * 0.5 * (1.0 - 2.0 * nu);
    }

    //	DB = D * B (3 x 6)
    double DB[3][6];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 6; j++)
        {
            DB[i][j] = 0.0;
            for (int k = 0; k < 3; k++)
                DB[i][j] += D[i][k] * B[k][j];
        }

    //	K = t * A * B^T * DB, only upper triangle
    double K[6][6] = { 0.0 };
    for (int i = 0; i < 6; i++)
        for (int j = i; j < 6; j++)
            for (int k = 0; k < 3; k++)
                K[i][j] += B[k][i] * DB[k][j];

    double coef = material->t * Area;

    int idx = 0;
    for (int j = 0; j < 6; j++)
        for (int i = j; i >= 0; i--)
            Matrix[idx++] = coef * K[i][j];
}

//	Calculate element stress (sx, sy, txy)
void CT3::ElementStress(double* stress, double* Displacement)
{
    //	Extract element displacements
    double u[6];
    for (int i = 0; i < 6; i++)
    {
        if (LocationMatrix_[i])
            u[i] = Displacement[LocationMatrix_[i] - 1];
        else
            u[i] = 0.0;
    }

    //	Node coordinates
    double x1 = nodes_[0]->XYZ[0];
    double y1 = nodes_[0]->XYZ[1];
    double x2 = nodes_[1]->XYZ[0];
    double y2 = nodes_[1]->XYZ[1];
    double x3 = nodes_[2]->XYZ[0];
    double y3 = nodes_[2]->XYZ[1];

    double b1 = y2 - y3;
    double b2 = y3 - y1;
    double b3 = y1 - y2;
    double c1 = x3 - x2;
    double c2 = x1 - x3;
    double c3 = x2 - x1;

    double Area2 = b1 * c2 - b2 * c1;
    double inv2A = 1.0 / Area2;

    //	B matrix
    double B[3][6];
    B[0][0] = b1 * inv2A;  B[0][1] = 0.0;        B[0][2] = b2 * inv2A;  B[0][3] = 0.0;        B[0][4] = b3 * inv2A;  B[0][5] = 0.0;
    B[1][0] = 0.0;         B[1][1] = c1 * inv2A;  B[1][2] = 0.0;         B[1][3] = c2 * inv2A;  B[1][4] = 0.0;         B[1][5] = c3 * inv2A;
    B[2][0] = c1 * inv2A;  B[2][1] = b1 * inv2A;  B[2][2] = c2 * inv2A;  B[2][3] = b2 * inv2A;  B[2][4] = c3 * inv2A;  B[2][5] = b3 * inv2A;

    //	Strain
    double e[3] = { 0.0 };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 6; j++)
            e[i] += B[i][j] * u[j];

    //	D matrix
    CT3Material* material = dynamic_cast <CT3Material*>(ElementMaterial_);
    double D[3][3];
    double E = material->E;
    double nu = material->nu;
    double factor;

    if (material->type == 1)	// Plane stress
    {
        factor = E / (1.0 - nu * nu);
        D[0][0] = factor;        D[0][1] = factor * nu;  D[0][2] = 0.0;
        D[1][0] = factor * nu;   D[1][1] = factor;       D[1][2] = 0.0;
        D[2][0] = 0.0;           D[2][1] = 0.0;          D[2][2] = factor * 0.5 * (1.0 - nu);
    }
    else	// Plane strain
    {
        factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        D[0][0] = factor * (1.0 - nu);  D[0][1] = factor * nu;          D[0][2] = 0.0;
        D[1][0] = factor * nu;          D[1][1] = factor * (1.0 - nu);  D[1][2] = 0.0;
        D[2][0] = 0.0;                  D[2][1] = 0.0;                  D[2][2] = factor * 0.5 * (1.0 - 2.0 * nu);
    }

    //	Stress = D * strain
    for (int i = 0; i < 3; i++)
    {
        stress[i] = 0.0;
        for (int j = 0; j < 3; j++)
            stress[i] += D[i][j] * e[j];
    }
}