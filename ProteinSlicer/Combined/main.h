/*
 *  main.h
 *  ProteinSlicer
 *
 *  Created by Frank Heinrich on 10/04/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#define MATRIXSIZE 140
#define UNITSIZE 0.5
#define PROBERADIUS 1.4 //1.4
#define RECURSIVECALLMAX 20000

#define NUMBEROFSTRUCTURES 1

// if BATCHMODE is ON (1) strPath and strFileNameRoot have to be used
// individual structure file names are expected to have the structure
// filenameroot_tiltxxx_oriexxx.pdb, with tilt being the beta Euler
// angle and orie being the gamma Euler angle for the rotated proteins.
// The output file will only differ in the ending .txt instead of .pdb

// if BATCHMODE is OFF (0) for single file processing, use filename and
// FilenameWrite
#define BATCHMODE 0
#define TILTSTART 0 //was 0
#define TILTEND 185  // was 185
#define TILTSTEP 5
#define ORIESTART 0  // was 0
#define ORIEEND 360  // was 360
#define ORIESTEP 5

#define COMPLEX 1
#define PROTCHAIN "A"
#define DEUTCHAIN "B"

char strPath[200]="/Users/frank/Documents/projects/proteins/KRas/StructureModeling/4G0N_RAS_RBD_complex/";
char strFileNameRoot[20]="protein";
char filename[200]="/Users/frank/Documents/projects/proteins/KRas/StructureModeling/KRas-RBD-Simanshu/rotated/bestfit_145tilt_40rot/protein_tilt145_orie40.pdb";
char FilenameWrite[200]="/Users/frank/Documents/projects/proteins/KRas/StructureModeling/KRas-RBD-Simanshu/rotated/bestfit_145tilt_40rot/protein_tilt145_orie40.txt";

// volume matrix starts from 0 positon to matrix -> ceil(pos/unit)
// value description:
//   0 - cell empty (void)
//   1 - cell with protein
//   2 - cell with solvent
//   4 - Help value for the recursive algorithm
//       in order to mark a bin, where it should have
//       stepped in with a fnFillProbeVolume() call but could
//       not, because of an already reached recursive depth limit.
//       The bin is filled with solvent.
int iMatrix[MATRIXSIZE][MATRIXSIZE][MATRIXSIZE];

// nsl matrix
// second argument is for protein without proton exchange (0)
// and protein with proton exchange (1)
double dMatrixnSL[MATRIXSIZE][2];

// storage matrix
// second argument is for protein without proton exchange (0)
// and protein with proton exchange (1)
// area [2]
double dMatrixStore[MATRIXSIZE][3];


double dPosXMin, dPosXMax, dPosYMin, dPosYMax, dPosZMin, dPosZMax;
double iRecursiveCallCounter;
long int iSolventAtomsFilled;
int iLastFraction, iProbeRadius, iProbeFrame;


enum { ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU,
    MET, ASN, PRO, GLN, ARG, SER, THR, VAL, TRP, TYR};

const char cNameAminoAcid[20][4]=
{"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"};

const double dnSLAminoAcid[20]=
{1.6344e-4, 1.9191e-4, 3.8343e-4, 3.7509e-4, 4.1268e-4,
    1.7178e-4, 4.7406e-4, 1.3842e-4, 1.5660e-4, 1.3842e-4,
    1.7523e-4, 3.4356e-4, 2.2158e-4, 3.3522e-4, 3.4260e-4,
    2.2149e-4, 2.1315e-4, 1.4676e-4, 6.0123e-4, 4.7073e-4};

const double dnSLDAminoAcid[20]=
{1.79e-6, 1.82e-6, 3.08e-6, 2.42e-6, 2.03e-6, 2.59e-6, 2.83e-6,
    8.20e-7, 9.14e-7, 8.24e-7, 1.03e-6, 2.54e-6, 1.71e-6, 2.08e-6,
    1.70e-6, 2.24e-6, 1.75e-6, 1.04e-6, 2.53e-6, 2.31e-6};

const double dVolumeAminoAcid[20]=
{91.5, 105.6, 124.5, 155.1, 203.4, 66.4, 167.3, 168.8, 171.3,
    167.9, 170.8, 135.2, 129.3, 161.1, 202.1, 99.1, 122.1, 141.7,
    237.6, 203.6};



const int iDimensionAtomPropertyTable = 28;
const double dAtomPropertyTable[28][3]=
{1, 5.805e-5, 1.60, //O doublebond C
    2, 5.805e-5, 1.70, //OH
    3, 5.805e-5, 1.60, //OCO
    4, 9.300e-5, 1.65, //NH
    5, 9.300e-5, 1.70, //NH2
    6, 9.300e-5, 1.75, //NH3
    7, 6.48e-5, 1.85, //CH
    8, 6.48e-5, 1.90, //CH2
    9, 6.48e-5, 1.95, //CH3
    10,	6.48e-5, 1.80, //CAr
    11, 6.48e-5, 1.90, //CHAR
    12, 2.847e-5, 1.90,	//S
    21, 3.630e-5, 1.50, //Na
    22, 4.700e-5, 1.50, //Ca
    23, 3.635e-5, 1.50, //Cr
    24, 5.070e-5, 1.50, //Ba
    25, 8.240e-5, 1.50, //La
    26, 7.630e-5, 1.50, //Au
    31, 6.48e-5, 2.00, //CAlNu
    32, 6.48e-5, 1.77, //CArNu
    33, 5.805e-5, 1.40, //OSuNu
    34, 5.805e-5, 1.64, //O doublebond CNu
    35, 5.805e-5, 1.64, //OPONu
    36, 9.300e-5, 1.55, //NArNu
    37, 9.300e-5, 1.86, //NAlNu
    38, 5.130e-5, 1.80, //P
    98, 6.674e-5, 1.20, //H exchangeable   must be second last element in table
    99, -3.741e-5, 1.20}; //H not exchangeable   must be last element in table

//profiling of all amino residues using Connolly's nomenclature
const int iDimensionAtomNameTable = 902;
const char cAtomNameTable[902][3][6]=
{
    "ALA", "N", "4",
    "ALA", "CA", "31",
    "ALA", "C", "31",
    "ALA", "O", "1",
    "ALA", "OXT", "2",    //O for terminal residue
    "ALA", "CB", "9",
    "ALA", "H", "98",
    "ALA", "HN", "98",    // terminal H
    "ALA", "HA", "99",
    "ALA", "1HB", "99",
    "ALA", "2HB", "99",
    "ALA", "3HB", "99",
    "ALA", "1H", "99",    //H for terminal residues
    "ALA", "2H", "99",
    "ALA", "3H", "99",
    "ALA", "HT1", "99",    //H for terminal residues
    "ALA", "HT2", "99",
    "ALA", "HT3", "99",
    "ARG", "N", "4",
    "ARG", "CA", "31",
    "ARG", "C", "31",
    "ARG", "O", "1",
    "ARG", "OXT", "2",
    "ARG", "CB", "8",
    "ARG", "CG", "8",
    "ARG", "CD", "8",
    "ARG", "NE", "4",
    "ARG", "CZ", "31",
    "ARG", "NH1", "37",
    "ARG", "NH2", "37",
    "ARG", "H", "98",
    "ARG", "HE", "98",
    "ARG", "1HH1", "98",
    "ARG", "2HH1", "98",
    "ARG", "1HH2", "98",
    "ARG", "2HH2", "98",
    "ARG", "HN", "98",    // terminal H
    "ARG", "HA", "99",
    "ARG", "1HB", "99",
    "ARG", "2HB", "99",
    "ARG", "3HB", "99",
    "ARG", "1HG", "99",
    "ARG", "2HG", "99",
    "ARG", "3HG", "99",
    "ARG", "1HD", "99",
    "ARG", "2HD", "99",
    "ARG", "3HD", "99",
    "ARG", "1H", "99",    //H for terminal residues
    "ARG", "2H", "99",
    "ARG", "3H", "99",
    "ARG", "HT1", "99",    //H for terminal residues
    "ARG", "HT2", "99",
    "ARG", "HT3", "99",
    "ASN", "N", "4",
    "ASN", "CA", "31",
    "ASN", "C", "31",
    "ASN", "O", "1",
    "ASN", "OXT", "2",
    "ASN", "OT1", "2",
    "ASN", "OT2", "2",
    "ASN", "CB", "8",
    "ASN", "CG", "31",
    "ASN", "OD1", "34",
    "ASN", "ND2", "31",
    "ASN", "H", "98",
    "ASN", "1HD2", "98",
    "ASN", "2HD2", "98",
    "ASN", "HN", "98",    // terminal H
    "ASN", "HA", "99",
    "ASN", "1HB", "99",
    "ASN", "2HB", "99",
    "ASN", "3HB", "99",
    "ASN", "1H", "99",    //H for terminal residues
    "ASN", "2H", "99",
    "ASN", "3H", "99",
    "ASN", "HT1", "99",    //H for terminal residues
    "ASN", "HT2", "99",
    "ASN", "HT3", "99",
    "ASP", "N", "4",
    "ASP", "CA", "31",
    "ASP", "C", "31",
    "ASP", "O", "1",
    "ASP", "OXT", "2",
    "ASP", "CB", "8",
    "ASP", "CG", "31",
    "ASP", "OD1", "34",
    "ASP", "OD1A", "34",
    "ASP", "OD1B", "34",
    "ASP", "OD2", "34",
    "ASP", "OD2A", "34",
    "ASP", "OD2B", "34",
    "ASP", "H", "98",
    "ASP", "HN", "98",    // terminal H
    "ASP", "HA", "99",
    "ASP", "1HB", "99",
    "ASP", "2HB", "99",
    "ASP", "HB2A", "99",
    "ASP", "HB2B", "99",
    "ASP", "3HB", "99",
    "ASP", "HB3A", "99",
    "ASP", "HB3B", "99",
    "ASP", "1H", "99",    //H for terminal residues
    "ASP", "2H", "99",
    "ASP", "3H", "99",
    "ASP", "HT1", "99",    //H for terminal residues
    "ASP", "HT2", "99",
    "ASP", "HT3", "99",
    "CYS", "N", "4",
    "CYS", "CA", "31",
    "CYS", "C", "31",
    "CYS", "O", "1",
    "CYS", "OXT", "2",
    "CYS", "1OT", "2",
    "CYS", "2OT", "2",
    "CYS", "CB", "8",
    "CYS", "SG", "12",
    "CYS", "H", "98",
    "CYS", "HG", "98",
    "CYS", "1HG", "98",
    "CYS", "HN", "98",    // terminal H
    "CYS", "HA", "99",
    "CYS", "1HB", "99",
    "CYS", "2HB", "99",
    "CYS", "3HB", "99",
    "CYS", "1H", "99",    //H for terminal residues
    "CYS", "2H", "99",
    "CYS", "3H", "99",
    "CYS", "HT1", "99",    //H for terminal residues
    "CYS", "HT2", "99",
    "CYS", "HT3", "99",
    "GLN", "N", "4",
    "GLN", "CA", "31",
    "GLN", "C", "31",
    "GLN", "O", "1",
    "GLN", "OXT", "2",
    "GLN", "OT1", "2",
    "GLN", "OT2", "2",
    "GLN", "CB", "8",
    "GLN", "CG", "8",
    "GLN", "CD", "31",
    "GLN", "OE1", "34",
    "GLN", "OE1A", "34",
    "GLN", "OE1B", "34",
    "GLN", "NE2", "31",
    "GLN", "NE2A", "31",
    "GLN", "NE2B", "31",
    "GLN", "H", "98",
    "GLN", "1HE2", "98",
    "GLN", "2HE2", "98",
    "GLN", "HE21A", "98",
    "GLN", "HE21B", "98",
    "GLN", "HE22A", "98",
    "GLN", "HE22B", "98",
    "GLN", "HN", "98",    // terminal H
    "GLN", "HA", "99",
    "GLN", "1HB", "99",
    "GLN", "2HB", "99",
    "GLN", "HB2A", "99",
    "GLN", "HB2B", "99",
    "GLN", "3HB", "99",
    "GLN", "HB3A", "99",
    "GLN", "HB3B", "99",
    "GLN", "1HG", "99",
    "GLN", "2HG", "99",
    "GLN", "HG2A", "99",
    "GLN", "HG2B", "99",
    "GLN", "3HG", "99",
    "GLN", "HG3A", "99",
    "GLN", "HG3B", "99",
    "GLN", "1H", "99",    //H for terminal residues
    "GLN", "2H", "99",
    "GLN", "3H", "99",
    "GLN", "HT1", "99",    //H for terminal residues
    "GLN", "HT2", "99",
    "GLN", "HT3", "99",
    "GLU", "N", "4",
    "GLU", "CA", "31",
    "GLU", "C", "31",
    "GLU", "O", "1",
    "GLU", "OXT", "2",
    "GLU", "CB", "8",
    "GLU", "CG", "8",
    "GLU", "CD", "31",
    "GLU", "OE1", "34",
    "GLU", "OE1A", "34",
    "GLU", "OE1B", "34",
    "GLU", "OE2", "34",
    "GLU", "OE2A", "34",
    "GLU", "OE2B", "34",
    "GLU", "H", "98",
    "GLU", "HN", "98",    // terminal H
    "GLU", "HA", "99",
    "GLU", "1HB", "99",
    "GLU", "2HB", "99",
    "GLU", "HB2A", "99",
    "GLU", "HB2B", "99",
    "GLU", "3HB", "99",
    "GLU", "HB3A", "99",
    "GLU", "HB3B", "99",
    "GLU", "1HG", "99",
    "GLU", "2HG", "99",
    "GLU", "HG2A", "99",
    "GLU", "HG2B", "99",
    "GLU", "3HG", "99",
    "GLU", "HG3A", "99",
    "GLU", "HG3B", "99",
    "GLU", "1H", "99",    //H for terminal residues
    "GLU", "2H", "99",
    "GLU", "3H", "99",
    "GLU", "HT1", "99",    //H for terminal residues
    "GLU", "HT2", "99",
    "GLU", "HT3", "99",
    "GLY", "N", "4",
    "GLY", "CA", "31",
    "GLY", "C", "31",
    "GLY", "O", "1",
    "GLY", "OXT", "2",
    "GLY", "CB", "9",
    "GLY", "H", "98",
    "GLY", "HN", "98",    // terminal H
    "GLY", "HA", "99",
    "GLY", "1HA", "99",
    "GLY", "2HA", "99",
    "GLY", "3HA", "99",
    "GLY", "1H", "99",    //H for terminal residues
    "GLY", "2H", "99",
    "GLY", "3H", "99",
    "GLY", "HT1", "99",    //H for terminal residues
    "GLY", "HT2", "99",
    "GLY", "HT3", "99",
    "HIS", "N", "4",
    "HIS", "CA", "31",
    "HIS", "C", "31",
    "HIS", "O", "1",
    "HIS", "OT1", "2",
    "HIS", "OT2", "2",
    "HIS", "OXT", "2",
    "HIS", "CB", "8",
    "HIS", "CG", "32",
    "HIS", "ND1", "4",
    "HIS", "CD2", "32",
    "HIS", "CE1", "32",
    "HIS", "NE2", "4",
    "HIS", "H", "98",
    "HIS", "HN", "98",    // terminal H
    "HIS", "HA", "99",
    "HIS", "1HB", "99",
    "HIS", "2HB", "99",
    "HIS", "3HB", "99",
    "HIS", "HD1", "98",
    "HIS", "HD2", "99",
    "HIS", "HE1", "99",
    "HIS", "HE2", "98",
    "HIS", "1H", "99",    //H for terminal residues
    "HIS", "2H", "99",
    "HIS", "3H", "99",
    "HIS", "HT1", "99",    //H for terminal residues
    "HIS", "HT2", "99",
    "HIS", "HT3", "99",
    "ILE", "N", "4",
    "ILE", "CA", "31",
    "ILE", "C", "31",
    "ILE", "O", "1",
    "ILE", "OXT", "2",
    "ILE", "CB", "7",
    "ILE", "CG1", "8",
    "ILE", "CG1A", "8",
    "ILE", "CG1B", "8",
    "ILE", "CG2", "9",
    "ILE", "CG2A", "8",
    "ILE", "CG2B", "8",
    "ILE", "CD1", "9",
    "ILE", "CD1A", "9",
    "ILE", "CD1B", "9",
    "ILE", "H", "98",
    "ILE", "HN", "98",    // terminal H
    "ILE", "HA", "99",
    "ILE", "HB", "99",
    "ILE", "1HG1", "99",
    "ILE", "2HG1", "99",
    "ILE", "HG12A", "99",
    "ILE", "HG12B", "99",
    "ILE", "1HG2", "99",
    "ILE", "HG21A", "99",
    "ILE", "HG21B", "99",
    "ILE", "2HG2", "99",
    "ILE", "HG22A", "99",
    "ILE", "HG22B", "99",
    "ILE", "HG22B", "99",
    "ILE", "3HG1", "99",
    "ILE", "HG13A", "99",
    "ILE", "HG13B", "99",
    "ILE", "3HG2", "99",
    "ILE", "HG23A", "99",
    "ILE", "HG23B", "99",
    "ILE", "1HD1", "99",
    "ILE", "HD11A", "99",
    "ILE", "HD11B", "99",
    "ILE", "2HD1", "99",
    "ILE", "HD12A", "99",
    "ILE", "HD12B", "99",
    "ILE", "3HD1", "99",
    "ILE", "HD13A", "99",
    "ILE", "HD13B", "99",
    "ILE", "1H", "99",    //H for terminal residues
    "ILE", "2H", "99",
    "ILE", "3H", "99",
    "ILE", "HT1", "99",    //H for terminal residues
    "ILE", "HT2", "99",
    "ILE", "HT3", "99",
    "LEU", "N", "4",
    "LEU", "CA", "31",
    "LEU", "C", "31",
    "LEU", "O", "1",
    "LEU", "OXT", "2",
    "LEU", "CB", "8",
    "LEU", "CG", "7",
    "LEU", "CD1", "9",
    "LEU", "CD2", "9",
    "LEU", "H", "98",
    "LEU", "HN", "98",    // terminal H
    "LEU", "HA", "99",
    "LEU", "1HB", "99",
    "LEU", "2HB", "99",
    "LEU", "3HB", "99",
    "LEU", "HG", "99",
    "LEU", "1HD1", "99",
    "LEU", "2HD1", "99",
    "LEU", "3HD1", "99",
    "LEU", "1HD2", "99",
    "LEU", "2HD2", "99",
    "LEU", "3HD2", "99",
    "LEU", "1H", "99",    //H for terminal residues
    "LEU", "2H", "99",
    "LEU", "3H", "99",
    "LEU", "HT1", "99",    //H for terminal residues
    "LEU", "HT2", "99",
    "LEU", "HT3", "99",
    "LYS", "N", "4",
    "LYS", "CA", "31",
    "LYS", "C", "31",
    "LYS", "O", "1",
    "LYS", "OXT", "2",
    "LYS", "OT1", "2",
    "LYS", "OT2", "2",
    "LYS", "CB", "8",
    "LYS", "CG", "8",
    "LYS", "CD", "8",
    "LYS", "CE", "31",
    "LYS", "NZ", "5",
    "LYS", "H", "98",
    "LYS", "1HZ", "98",
    "LYS", "2HZ", "98",
    "LYS", "3HZ", "98",
    "LYS", "HN", "98",    // terminal H
    "LYS", "HA", "99",
    "LYS", "1HB", "99",
    "LYS", "2HB", "99",
    "LYS", "3HB", "99",
    "LYS", "1HG", "99",
    "LYS", "2HG", "99",
    "LYS", "3HG", "99",
    "LYS", "1HD", "99",
    "LYS", "2HD", "99",
    "LYS", "3HD", "99",
    "LYS", "1HE", "99",
    "LYS", "2HE", "99",
    "LYS", "3HE", "99",
    "LYS", "1H", "99",    //H for terminal residues
    "LYS", "2H", "99",
    "LYS", "3H", "99",
    "LYS", "HT1", "99",    //H for terminal residues
    "LYS", "HT2", "99",
    "LYS", "HT3", "99",
    "MET", "N", "4",
    "MET", "CA", "31",
    "MET", "C", "31",
    "MET", "O", "1",
    "MET", "OXT", "2",
    "MET", "CB", "8",
    "MET", "CG", "8",
    "MET", "SD", "12",
    "MET", "CE", "9",
    "MET", "H", "98",
    "MET", "HN", "98",    // terminal H
    "MET", "HA", "99",
    "MET", "1HB", "99",
    "MET", "2HB", "99",
    "MET", "3HB", "99",
    "MET", "1HG", "99",
    "MET", "2HG", "99",
    "MET", "3HG", "99",
    "MET", "1HE", "99",
    "MET", "2HE", "99",
    "MET", "3HE", "99",
    "MET", "1H", "99",    //H for terminal residues
    "MET", "2H", "99",
    "MET", "3H", "99",
    "MET", "HT1", "99",    //H for terminal residues
    "MET", "HT2", "99",
    "MET", "HT3", "99",
    "MSE", "N", "4",
    "MSE", "CA", "31",
    "MSE", "C", "31",
    "MSE", "O", "1",
    "MSE", "OXT", "2",
    "MSE", "CB", "8",
    "MSE", "CG", "8",
    "MSE", "SE", "12", //treat the SE as sulphur, what it actually is
    "MSE", "CE", "9",
    "MSE", "H", "98",
    "MSE", "HN", "98",    // terminal H
    "MSE", "HA", "99",
    "MSE", "1HB", "99",
    "MSE", "2HB", "99",
    "MSE", "3HB", "99",
    "MSE", "1HG", "99",
    "MSE", "2HG", "99",
    "MSE", "3HG", "99",
    "MSE", "1HE", "99",
    "MSE", "2HE", "99",
    "MSE", "3HE", "99",
    "MSE", "HB1", "99",
    "MSE", "HB2", "99",
    "MSE", "HB3", "99",
    "MSE", "HG1", "99",
    "MSE", "HG2", "99",
    "MSE", "HG3", "99",
    "MSE", "HE1", "99",
    "MSE", "HE2", "99",
    "MSE", "HE3", "99",
    "MSE", "1H", "99",    //H for terminal residues
    "MSE", "2H", "99",
    "MSE", "3H", "99",
    "MSE", "HT1", "99",    //H for terminal residues
    "MSE", "HT2", "99",
    "MSE", "HT3", "99",
    "MYR", "O1", "1",
    "MYR", "C1", "9",
    "MYR", "C2", "8",
    "MYR", "C3", "8",
    "MYR", "C4", "8",
    "MYR", "C5", "8",
    "MYR", "C6", "8",
    "MYR", "C7", "8",
    "MYR", "C8", "8",
    "MYR", "C9", "8",
    "MYR", "C10", "8",
    "MYR", "C11", "8",
    "MYR", "C12", "8",
    "MYR", "C13", "8",
    "MYR", "C14", "8",
    "MYR", "H2A", "99",
    "MYR", "H2B", "99",
    "MYR", "1H2", "99",
    "MYR", "2H2", "99",
    "MYR", "H3A", "99",
    "MYR", "H3B", "99",
    "MYR", "1H3", "99",
    "MYR", "2H3", "99",
    "MYR", "H4A", "99",
    "MYR", "H4B", "99",
    "MYR", "1H4", "99",
    "MYR", "2H4", "99",
    "MYR", "H5A", "99",
    "MYR", "H5B", "99",
    "MYR", "1H5", "99",
    "MYR", "2H5", "99",
    "MYR", "H6A", "99",
    "MYR", "H6B", "99",
    "MYR", "1H6", "99",
    "MYR", "2H6", "99",
    "MYR", "H7A", "99",
    "MYR", "H7B", "99",
    "MYR", "1H7", "99",
    "MYR", "2H7", "99",
    "MYR", "H8A", "99",
    "MYR", "H8B", "99",
    "MYR", "1H8", "99",
    "MYR", "2H8", "99",
    "MYR", "H9A", "99",
    "MYR", "H9B", "99",
    "MYR", "1H9", "99",
    "MYR", "2H9", "99",
    "MYR", "H10A", "99",
    "MYR", "H10B", "99",
    "MYR", "1H10", "99",
    "MYR", "2H10", "99",
    "MYR", "H11A", "99",
    "MYR", "H11B", "99",
    "MYR", "1H11", "99",
    "MYR", "2H11", "99",
    "MYR", "H12A", "99",
    "MYR", "H12B", "99",
    "MYR", "1H12", "99",
    "MYR", "2H12", "99",
    "MYR", "H13A", "99",
    "MYR", "H13B", "99",
    "MYR", "1H13", "99",
    "MYR", "2H13", "99",
    "MYR", "H14A", "99",
    "MYR", "H14B", "99",
    "MYR", "H14C", "99",
    "MYR", "1H14", "99",
    "MYR", "2H14", "99",
    "MYR", "H101", "99",
    "MYR", "H102", "99",
    "MYR", "H121", "99",
    "MYR", "H122", "99",
    "MYR", "H131", "99",
    "MYR", "H132", "99",
    "MYR", "H141", "99",
    "MYR", "H142", "99",
    "MYR", "H143", "99",
    "NON", "C", "8",      // NON = no residue given in file
    "NON", "O", "2",      // atom types chosen to match the most abundant
    "NON", "S", "12",
    "NON", "H", "99",
    "NON", "D", "98",
    "NON", "P", "38",
    "NON", "N", "4",
    "PHE", "N", "4",
    "PHE", "CA", "31",
    "PHE", "C", "31",
    "PHE", "O", "1",
    "PHE", "OXT", "2",
    "PHE", "CB", "8",
    "PHE", "CG", "10",
    "PHE", "CD1", "11",
    "PHE", "CD2", "11",
    "PHE", "CE1", "11",
    "PHE", "CE2", "11",
    "PHE", "CZ", "11",
    "PHE", "H", "98",
    "PHE", "HN", "98",    // terminal H
    "PHE", "HA", "99",
    "PHE", "1HB", "99",
    "PHE", "2HB", "99",
    "PHE", "3HB", "99",
    "PHE", "HD1", "99",
    "PHE", "HD2", "99",
    "PHE", "HE1", "99",
    "PHE", "HE2", "99",
    "PHE", "HZ", "99",
    "PHE", "1H", "99",    //H for terminal residues
    "PHE", "2H", "99",
    "PHE", "3H", "99",
    "PHE", "HT1", "99",    //H for terminal residues
    "PHE", "HT2", "99",
    "PHE", "HT3", "99",
    "PRO", "N", "4",
    "PRO", "CA", "31",
    "PRO", "C", "31",
    "PRO", "O", "1",
    "PRO", "OXT", "2",
    "PRO", "CB", "8",
    "PRO", "CG", "8",
    "PRO", "CD", "31",
    "PRO", "HN", "98",    // terminal H
    "PRO", "H", "99",
    "PRO", "HA", "99",
    "PRO", "1HB", "99",
    "PRO", "2HB", "99",
    "PRO", "3HB", "99",
    "PRO", "1HG", "99",
    "PRO", "2HG", "99",
    "PRO", "3HG", "99",
    "PRO", "1HD", "99",
    "PRO", "2HD", "99",
    "PRO", "3HD", "99",
    "PRO", "1H", "99",    //H for terminal residues
    "PRO", "2H", "99",
    "PRO", "3H", "99",
    "PRO", "HT1", "99",    //H for terminal residues
    "PRO", "HT2", "99",
    "PRO", "HT3", "99",
    "SER", "N", "4",
    "SER", "CA", "31",
    "SER", "C", "31",
    "SER", "O", "1",
    "SER", "OXT", "2",
    "SER", "CB", "31",
    "SER", "OG", "2",
    "SER", "H", "98",
    "SER", "HG", "98",
    "SER", "1HG", "98",
    "SER", "HN", "98",    // terminal H
    "SER", "HA", "99",
    "SER", "1HB", "99",
    "SER", "2HB", "99",
    "SER", "HB2A", "99",
    "SER", "3HB", "99",
    "SER", "HB3A", "99",
    "SER", "1H", "99",    //H for terminal residues
    "SER", "2H", "99",
    "SER", "3H", "99",
    "SER", "HT1", "99",    //H for terminal residues
    "SER", "HT2", "99",
    "SER", "HT3", "99",
    "THR", "N", "4",
    "THR", "CA", "31",
    "THR", "C", "31",
    "THR", "O", "1",
    "THR", "OXT", "2",
    "THR", "CB", "31",
    "THR", "OG1", "2",
    "THR", "CG2", "9",
    "THR", "H", "98",
    "THR", "HG1", "98",
    "THR", "HN", "98",    // terminal H
    "THR", "HA", "99",
    "THR", "HB", "99",
    "THR", "1HG2", "99",
    "THR", "2HG2", "99",
    "THR", "3HG2", "99",
    "THR", "1H", "99",    //H for terminal residues
    "THR", "2H", "99",
    "THR", "3H", "99",
    "THR", "HT1", "99",    //H for terminal residues
    "THR", "HT2", "99",
    "THR", "HT3", "99",
    "TRP", "N", "4",
    "TRP", "CA", "31",
    "TRP", "C", "31",
    "TRP", "O", "1",
    "TRP", "OXT", "2",
    "TRP", "CB", "8",
    "TRP", "CG", "10",
    "TRP", "CD1", "32",
    "TRP", "CD2", "10",
    "TRP", "NE1", "4",
    "TRP", "CE2", "32",
    "TRP", "CE3", "11",
    "TRP", "CZ2", "11",
    "TRP", "CZ3", "11",
    "TRP", "H", "98",
    "TRP", "HE1", "98",
    "TRP", "HN", "98",    // terminal H
    "TRP", "HA", "99",
    "TRP", "1HB", "99",
    "TRP", "2HB", "99",
    "TRP", "3HB", "99",
    "TRP", "HD1", "99",
    "TRP", "HE3", "99",
    "TRP", "HZ2", "99",
    "TRP", "HZ3", "99",
    "TRP", "HH2", "99",
    "TRP", "CH2", "11",
    "TRP", "1H", "99",    //H for terminal residues
    "TRP", "2H", "99",
    "TRP", "3H", "99",
    "TRP", "HT1", "99",    //H for terminal residues
    "TRP", "HT2", "99", 
    "TRP", "HT3", "99",
    "TYR", "N", "4",
    "TYR", "NT", "4",
    "TYR", "CA", "31",
    "TYR", "CAT", "31",
    "TYR", "C", "31",
    "TYR", "O", "1",
    "TYR", "OXT", "2",
    "TYR", "CB", "8",
    "TYR", "CG", "10",
    "TYR", "CD1", "11",
    "TYR", "CD2", "11",
    "TYR", "CE1", "11",
    "TYR", "CE2", "11",
    "TYR", "CZ", "32",
    "TYR", "H", "98",
    "TYR", "HH", "98",
    "TYR", "OH", "2",
    "TYR", "HN", "98",    // terminal H
    "TYR", "HNT", "98",    // terminal H
    "TYR", "HA", "99",
    "TYR", "1HB", "99",
    "TYR", "2HB", "99",
    "TYR", "3HB", "99",
    "TYR", "HD1", "99",
    "TYR", "HD2", "99",
    "TYR", "HE1", "99",
    "TYR", "HE2", "99",
    "TYR", "1H", "99",    //H for terminal residues
    "TYR", "2H", "99", 
    "TYR", "3H", "99",
    "TYR", "HT1", "99",    //H for terminal residues
    "TYR", "HT2", "99", 
    "TYR", "HT3", "99",
    "VAL", "N", "4",
    "VAL", "CA", "31",
    "VAL", "C", "31",
    "VAL", "O", "1",
    "VAL", "OXT", "2",
    "VAL", "CB", "7",
    "VAL", "CG1", "9",
    "VAL", "CG1A", "9",
    "VAL", "CG1B", "9",
    "VAL", "CG2", "9",
    "VAL", "CG2A", "9",
    "VAL", "CG2B", "9",
    "VAL", "H", "98",
    "VAL", "HN", "98",    // terminal H
    "VAL", "HA", "99",
    "VAL", "HB", "99",
    "VAL", "1HG1", "99",
    "VAL", "HG11A", "99",
    "VAL", "HG11B", "99",
    "VAL", "2HG1", "99",
    "VAL", "HG12A", "99",
    "VAL", "HG12B", "99",
    "VAL", "3HG1", "99",
    "VAL", "HG13A", "99",
    "VAL", "HG13B", "99",
    "VAL", "1HG2", "99",
    "VAL", "HG21A", "99",
    "VAL", "HG21B", "99",
    "VAL", "2HG2", "99",
    "VAL", "HG22A", "99",
    "VAL", "HG22B", "99",
    "VAL", "3HG2", "99",
    "VAL", "HG23A", "99",
    "VAL", "HG23B", "99",
    "VAL", "1H", "99",    //H for terminal residues
    "VAL", "2H", "99", 
    "VAL", "3H", "99",
    "VAL", "HT1", "99",    //H for terminal residues
    "VAL", "HT2", "99", 
    "VAL", "HT3", "99",
    "AMA", "C1", "8",      //Sugar groups
    "AMA", "C2", "8",
    "AMA", "C3", "8",
    "AMA", "C4", "8",
    "AMA", "C5", "8",
    "AMA", "C6", "8",
    "AMA", "O1", "2",
    "AMA", "O2", "2",
    "AMA", "O3", "2",
    "AMA", "O4", "2",
    "AMA", "O5", "2",
    "AMA", "O6", "2",
    "AMA", "H1", "99",
    "AMA", "H2", "99",
    "AMA", "H3", "99",
    "AMA", "H4", "99",
    "AMA", "H5", "99",
    "AMA", "HO2", "98",
    "AMA", "HO3", "98",
    "AMA", "HO4", "98",
    "AMA", "HO6", "98",
    "AMA", "H61", "99",
    "AMA", "H62", "99",
    "BMA", "C1", "8",      //Sugar groups
    "BMA", "C2", "8",
    "BMA", "C3", "8",
    "BMA", "C4", "8",
    "BMA", "C5", "8",
    "BMA", "C6", "8",
    "BMA", "O1", "2",
    "BMA", "O2", "2",
    "BMA", "O3", "2",
    "BMA", "O4", "2",
    "BMA", "O5", "2",
    "BMA", "O6", "2",
    "BMA", "H1", "99",
    "BMA", "H2", "99",
    "BMA", "H3", "99",
    "BMA", "H4", "99",
    "BMA", "H5", "99",
    "BMA", "HO2", "98",
    "BMA", "HO3", "98",
    "BMA", "HO4", "98",
    "BMA", "HO6", "98",
    "BMA", "H61", "99",
    "BMA", "H62", "99",
    "BFU", "C1", "8",      //Sugar groups
    "BFU", "C2", "8",
    "BFU", "C3", "8",
    "BFU", "C4", "8",
    "BFU", "C5", "8",
    "BFU", "C6", "8",
    "BFU", "O1", "2",
    "BFU", "O2", "2",
    "BFU", "O3", "2",
    "BFU", "O4", "2",
    "BFU", "O5", "2",
    "BFU", "O6", "2",
    "BFU", "H1", "99",
    "BFU", "H2", "99",
    "BFU", "H3", "99",
    "BFU", "H4", "99",
    "BFU", "H5", "99",
    "BFU", "HO2", "98",
    "BFU", "HO3", "98",
    "BFU", "HO4", "98",
    "BFU", "HO6", "98",
    "BFU", "H61", "99",
    "BFU", "H62", "99",
    "BFU", "H63", "99",
    "BNA", "C1", "8",      //Sugar groups
    "BNA", "C2", "8",
    "BNA", "C3", "8",
    "BNA", "C4", "8",
    "BNA", "C5", "8",
    "BNA", "C6", "8",
    "BNA", "CT", "9",
    "BNA", "O1", "2",
    "BNA", "O2", "2",
    "BNA", "O3", "2",
    "BNA", "O4", "2",
    "BNA", "O5", "2",
    "BNA", "O6", "2",
    "BNA", "H1", "99",
    "BNA", "H2", "99",
    "BNA", "H3", "99",
    "BNA", "H4", "99",
    "BNA", "H5", "99",
    "BNA", "HO2", "98",
    "BNA", "HO3", "98",
    "BNA", "HO4", "98",
    "BNA", "HO6", "98",
    "BNA", "H51", "99",
    "BNA", "H52", "99",
    "BNA", "H61", "99",
    "BNA", "H62", "99",
    "BNA", "HT1", "99",
    "BNA", "HT2", "99",
    "BNA", "HT3", "99",
    "BNA", "N", "4",
    "BNA", "HN", "99",
    "BXY", "C1", "8",      //Sugar groups
    "BXY", "C2", "8",
    "BXY", "C3", "8",
    "BXY", "C4", "8",
    "BXY", "C5", "8",
    "BXY", "O1", "2",
    "BXY", "O2", "2",
    "BXY", "O3", "2",
    "BXY", "O4", "2",
    "BXY", "O5", "2",
    "BXY", "H1", "99",
    "BXY", "H2", "99",
    "BXY", "H3", "99",
    "BXY", "H4", "99",
    "BXY", "HO2", "98",
    "BXY", "HO3", "98",
    "BXY", "HO4", "98",
    "BXY", "H51", "99",
    "BXY", "H52", "99",
    "MAN", "C1", "8",      //Sugar groups
    "MAN", "C2", "8",
    "MAN", "C3", "8",
    "MAN", "C4", "8",
    "MAN", "C5", "8",
    "MAN", "C6", "8",
    "MAN", "O1", "2",
    "MAN", "O2", "2",
    "MAN", "O3", "2",
    "MAN", "O4", "2",
    "MAN", "O5", "2",
    "MAN", "O6", "2",
    "MAN", "H1", "99",
    "MAN", "H2", "99",
    "MAN", "H3", "99",
    "MAN", "H4", "99",
    "MAN", "H5", "99",
    "MAN", "HO2", "98",
    "MAN", "HO3", "98",
    "MAN", "HO4", "98",
    "MAN", "HO6", "98",
    "MAN", "H61", "99",
    "MAN", "H62", "99",
    "NAG", "C1", "8",      //Sugar groups
    "NAG", "C2", "8",
    "NAG", "C3", "8",
    "NAG", "C4", "8",
    "NAG", "C5", "8",
    "NAG", "C6", "8",
    "NAG", "C7", "8",
    "NAG", "C8", "8",
    "NAG", "N2", "4",
    "NAG", "O1", "2",
    "NAG", "O2", "2",
    "NAG", "O3", "2",
    "NAG", "O4", "2",
    "NAG", "O5", "2",
    "NAG", "O6", "2",
    "NAG", "O7", "2",
    "NAG", "H1", "99",
    "NAG", "H2", "99",
    "NAG", "H3", "99",
    "NAG", "H4", "99",
    "NAG", "H5", "99",
    "NAG", "HO", "98",
    "NAG", "HN", "98",
    "NAG", "HO2", "98",
    "NAG", "HO3", "98",
    "NAG", "HO4", "98",
    "NAG", "HO6", "98",
    "NAG", "H61", "99",
    "NAG", "H62", "99",
    "NAG", "H81", "99",
    "NAG", "H82", "99",
    "NAG", "H83", "99",
};


void fnAnalyzeMatrix();
double fnArrayPosDifferenz(int iX,int iY,int iZ,double dPosX,double dPosY,double dPosZ);
void fnCoord2Array(int* iX, int* iY, int* iZ, double dPosX, double dPosY, double dPosZ);
void fnCoord2ArraynSL(int* iZ, double dPosZ);

void fnFillnSlProfile(double dPosZ, char *cResidue, char *cAtomName, char *cChainIdentifier);           //filling nSL profile
void fnFillVanDerWaals(double dPosX,double dPosY,double iPosZ,char *cAtomName, char *cResidue);

int  fnDetermineLimits(char* strFilename);
void fnLoadFile(char* strFilename, int structure);

int  strcpprot(char *cAtomName, char *cTableEntry);

void fnFillSolvent(int iStartX, int iStartY, int iStartZ);
void fnFillProbeVolume(int iX, int iY, int iZ);

void fnStoreResults(int structure);

void fnWriteVolumeSliced(char* strFilename);
