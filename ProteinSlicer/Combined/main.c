#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "main.h"

//-------------------------------------------------------------------------------------------------------------------------
void fnAnalyzeMatrix()			//calculates statistic from matrix
{
    int iProteinCounter, iSolventCounter, iVoidCounter;
    int iX, iY, iZ, iZnSL;
    double fProtnSLSum, fDeutnSLSum;
    
    iProteinCounter=0; iSolventCounter=0; iVoidCounter=0;   //print intermediate result
    iZnSL=0; fProtnSLSum=0; fDeutnSLSum=0;
    for(iZ=0; iZ<MATRIXSIZE; iZ++)
    {
        for(iY=0; iY<MATRIXSIZE; iY++)
            for(iX=0; iX<MATRIXSIZE; iX++)
            {
                switch(iMatrix[iX][iY][iZ])
                {
                    case 0: iVoidCounter++; break;
                    case 1: iProteinCounter++; break;
                    case 2: iSolventCounter++; break;
                }
            }
        fProtnSLSum+=dMatrixnSL[iZnSL][0];
        fDeutnSLSum+=dMatrixnSL[iZnSL][1];
        iZnSL+=1;
    }
    
    printf("Void    : %lf Percent \n", (double)(iVoidCounter)/pow(MATRIXSIZE,3)*100);
    printf("Protein : %lf Percent \n", (double)(iProteinCounter)/pow(MATRIXSIZE,3)*100);
    printf("Solvent : %lf Percent \n", (double)(iSolventCounter)/pow(MATRIXSIZE,3)*100);
    printf("Total protein volume [Å^3]: %e \n", (double)(iProteinCounter+iVoidCounter)*pow(UNITSIZE,3));
    printf("prot nSL [Å^1]: %lf  \n", fProtnSLSum);
    printf("deut nSL [Å^1]: %lf  \n \n", fDeutnSLSum);
    printf("prot nSLD [10^-6 Å^-2]: %lf  \n", fProtnSLSum/(double)(iProteinCounter+iVoidCounter)/pow(UNITSIZE,3)*pow(10,6));
    printf("deut nSLD [10^-6 Å^-2]: %lf  \n \n", fDeutnSLSum/(double)(iProteinCounter+iVoidCounter)/pow(UNITSIZE,3)*pow(10,6));
    
    
}
//-------------------------------------------------------------------------------------------------------------------------


double fnArrayPosDifferenz(int iX,int iY,int iZ, double dPosX, double dPosY,double dPosZ)
{  //calculate the distance between the middle of an array element and a 
    //position given in coordinates
    double dx, dy, dz;
    
    dx=(((double)(iX-iProbeFrame))*UNITSIZE+dPosXMin-dPosX);
    dy=(((double)(iY-iProbeFrame))*UNITSIZE+dPosYMin-dPosY);
    dz=(((double)(iZ-iProbeFrame))*UNITSIZE+dPosZMin-dPosZ);
    
    return sqrt(dx*dx+dy*dy+dz*dz);
    
}

//-------------------------------------------------------------------------------------------------------------------------

void fnCoord2ArraynSL(int* iZ, double dPosZ)
{	//converting a position given in coordinates to array indizees
	//using the minimum Position dPosZMin
    
    *iZ = (int)((dPosZ-dPosZMin)/UNITSIZE+0.5)+iProbeFrame;
}

void fnCoord2Array(int* iX, int* iY, int* iZ, double dPosX, double dPosY, double dPosZ)
{	//converting a position given in coordinates to array indizees
	//using the minimum Position XMin, YMin and ZMin
    
    *iX = (int)((dPosX-dPosXMin)/UNITSIZE+0.5)+iProbeFrame;
    *iY = (int)((dPosY-dPosYMin)/UNITSIZE+0.5)+iProbeFrame;
    *iZ = (int)((dPosZ-dPosZMin)/UNITSIZE+0.5)+iProbeFrame;
}

//-------------------------------------------------------------------------------------------------------------------------


void fnFillVanDerWaals(double dPosX,double dPosY,double dPosZ,char *cAtomName, char *cResidue)
{
    double dVanDerWaalsRadiusPlus;
    int iStartX, iStartY, iStartZ, iEndX, iEndY, iEndZ, iX, iY, iZ, iFoundInAtomNameTable, iFoundInAtomPropertyTable, i, i2, i3;
    
    //averaged Connolly's van der Waals radi
    //this script discriminates only between atom types
    iFoundInAtomNameTable=0;
    for (i=0; i<iDimensionAtomNameTable; i++)   //loop over whole atom/resiude list
    {
        if ((strcmp(cResidue,cAtomNameTable[i][0]) == 0) && (strcpprot(cAtomName,cAtomNameTable[i][1])==0))   //look for match in table
        {
            iFoundInAtomNameTable=1;
            iFoundInAtomPropertyTable=0;				
            i2=atoi(cAtomNameTable[i][2]);
            for (i3=0; i3<iDimensionAtomPropertyTable; i3++)   //find element in property table
                if (i2==dAtomPropertyTable[i3][0]) 
                {
                    iFoundInAtomPropertyTable=1;
                    dVanDerWaalsRadiusPlus=dAtomPropertyTable[i3][2];			//look up Connolly's van der Waals Radius
                    
                    i=iDimensionAtomNameTable;					//end loop
                }
            if (iFoundInAtomPropertyTable==0) {printf("Unknown atom indentifier %i \n", i2);}
        }
    }
    if (iFoundInAtomNameTable==0) 
    {
        printf("Unknown Atom %s %s \n", cResidue, cAtomName);
    }
    else {
    //    printf("Known Atom %s %s \n", cResidue, cAtomName);
    
    // calculates cube in which the van der Waals sphere is embedded
    fnCoord2Array(&iStartX, &iStartY, &iStartZ, dPosX-dVanDerWaalsRadiusPlus, dPosY-dVanDerWaalsRadiusPlus, dPosZ-dVanDerWaalsRadiusPlus);
    fnCoord2Array(&iEndX, &iEndY, &iEndZ, dPosX+dVanDerWaalsRadiusPlus, dPosY+dVanDerWaalsRadiusPlus, dPosZ+dVanDerWaalsRadiusPlus);
    
    //check if atom surrounding lies within array
    if (iStartX<0) {iStartX=0;}				
    if (iStartY<0) {iStartY=0;}
    if (iStartZ<0) {iStartZ=0;}
    if (iEndX>=MATRIXSIZE) {iEndX=MATRIXSIZE-1;}
    if (iEndY>=MATRIXSIZE) {iEndY=MATRIXSIZE-1;}
    if (iEndZ>=MATRIXSIZE) {iEndZ=MATRIXSIZE-1;}
    
    //mark all occupied cells in matrix
    for (iX=iStartX; iX<=iEndX; iX++)
        for (iY=iStartY; iY<=iEndY; iY++)
            for (iZ=iStartZ; iZ<=iEndZ; iZ++)
                if (fnArrayPosDifferenz(iX,iY,iZ,dPosX,dPosY,dPosZ)<=dVanDerWaalsRadiusPlus) 
                {
                    iMatrix[iX][iY][iZ]=1;
                }
    }
}

//-------------------------------------------------------------------------------------------------------------------------


void fnFillProbeVolume(int iX, int iY, int iZ)
{
    int iStartX, iStartY, iStartZ, iEndX, iEndY, iEndZ;
    int iXl, iYl, iZl, iProbeCollision, iVoidExists;
    
    iRecursiveCallCounter++;
    
    iStartX=iX-iProbeRadius; //limits of cube surrounding the
    iStartY=iY-iProbeRadius; //probe sphere
    iStartZ=iZ-iProbeRadius;
    iEndX=iX+iProbeRadius;
    iEndY=iY+iProbeRadius;
    iEndZ=iZ+iProbeRadius;
	
    if (iStartX<0) iStartX=0;   //check if limits are within array 
    if (iStartY<0) iStartY=0;
    if (iStartZ<0) iStartZ=0;
    if (iEndX>=MATRIXSIZE-1) iEndX=MATRIXSIZE-1;
    if (iEndY>=MATRIXSIZE-1) iEndY=MATRIXSIZE-1;
    if (iEndZ>=MATRIXSIZE-1) iEndZ=MATRIXSIZE-1;
	
    iProbeCollision=0; iVoidExists=0;				//check if sphere around start point is not occupied by protein
    for (iXl=iStartX; iXl<=iEndX; iXl++)			//or already completely filled with solvent
        for (iYl=iStartY; iYl<=iEndY; iYl++)
            for (iZl=iStartZ; iZl<=iEndZ; iZl++)		//check if bin is within probe sphere
                if (((double)((iX-iXl)*(iX-iXl)+(iY-iYl)*(iY-iYl)+(iZ-iZl)*(iZ-iZl)))<=((PROBERADIUS/UNITSIZE)*(PROBERADIUS/UNITSIZE))) 
                {
                    if (iMatrix[iXl][iYl][iZl]==0) {iVoidExists=1;}
                    if (iMatrix[iXl][iYl][iZl]==1) {iProbeCollision=1; iXl=iEndX; iYl=iEndY; iZl=iEndZ;}
                }
	
    
    if ((iProbeCollision==0) && (iVoidExists==1))			//if no protein collision and still empty bins, then fill volume
    {
        for (iXl=iStartX; iXl<=iEndX; iXl++)
            for (iYl=iStartY; iYl<=iEndY; iYl++)
                for (iZl=iStartZ; iZl<=iEndZ; iZl++)
                    if (((double)((iX-iXl)*(iX-iXl)+(iY-iYl)*(iY-iYl)+(iZ-iZl)*(iZ-iZl)))<=((PROBERADIUS/UNITSIZE)*(PROBERADIUS/UNITSIZE)))  
                    {
                        if (iMatrix[iXl][iYl][iZl]==0)						//cell empty?
                        {
                            iMatrix[iXl][iYl][iZl]=2;							//fill with solvent
                            iSolventAtomsFilled++;							//monitor overall progress for user
                            if (iSolventAtomsFilled/10000>iLastFraction)		//print statement every 10000 filled atoms
                            {
                                //printf("Total volume filled: %lf Percent \n", 
                                //       (double)(iSolventAtomsFilled)/(MATRIXSIZE*MATRIXSIZE*MATRIXSIZE)*100);
                                iLastFraction++;
                            } 
                        }
                    }
        
        
        
        //recursive calls into all six directions of space
        //simulating the moving of the probe sphere
        if ((iX-1)>=0) 
        {
            if (iRecursiveCallCounter<RECURSIVECALLMAX) {fnFillProbeVolume(iX-1,iY,iZ);}		//recursive depth limit reached?
            else {iMatrix[iX-1][iY][iZ]=4;}													//mark cell if recursive jump not possible
        }
        if ((iY-1)>=0) 
        {
            if (iRecursiveCallCounter<RECURSIVECALLMAX) {fnFillProbeVolume(iX,iY-1,iZ);}
            else {iMatrix[iX][iY-1][iZ]=4;}
        }
        if ((iZ-1)>=0) 
        {
            if (iRecursiveCallCounter<RECURSIVECALLMAX) {fnFillProbeVolume(iX,iY,iZ-1);}
            else {iMatrix[iX][iY][iZ-1]=4;}
        }
        if ((iX+1)<MATRIXSIZE) 
        {
            if (iRecursiveCallCounter<RECURSIVECALLMAX) {fnFillProbeVolume(iX+1,iY,iZ);}
            else {iMatrix[iX+1][iY][iZ]=4;}
        }
        if ((iY+1)<MATRIXSIZE) 
        {
            if (iRecursiveCallCounter<RECURSIVECALLMAX) {fnFillProbeVolume(iX,iY+1,iZ);}
            else {iMatrix[iX][iY+1][iZ]=4;}
        }
        if ((iZ+1)<MATRIXSIZE) 
        {
            if (iRecursiveCallCounter<RECURSIVECALLMAX) {fnFillProbeVolume(iX,iY,iZ+1);}
            else {iMatrix[iX][iY][iZ+1]=4;}
        }
    }
    
    iRecursiveCallCounter--;
}

//-------------------------------------------------------------------------------------------------------------------------

void fnFillnSlProfile(double dPosZ, char *cResidue, char *cAtomName, char *cChainIdentifier)  //filling nSL profile
{  
    double dnSLprot, dnSLdeut;
    int iPosZ,i, i2, i3, iFoundInAtomNameTable, iFoundInAtomPropertyTable;
    
    //printf("%s %lf \n", cAtomType, dPosZ);
    //printf("%d ",iNumber);
    
    //nSL from Mathias' table
    
    iFoundInAtomNameTable=0;
    for (i=0; i<iDimensionAtomNameTable; i++)   //loop over whole atom/resiude list
    {
        if ((strcmp(cResidue,cAtomNameTable[i][0]) == 0) && (strcpprot(cAtomName,cAtomNameTable[i][1])== 0))   //look for match in table
        {
            iFoundInAtomNameTable=1;
            if (strcmp(cAtomNameTable[i][2],"98") == 0)								//Exchangeable Hydrogen! Calculating nSL from deuteration ratio
            {
                dnSLprot=dAtomPropertyTable[iDimensionAtomPropertyTable-1][1];
                dnSLdeut=dAtomPropertyTable[iDimensionAtomPropertyTable-2][1];
                //printf("nSLFill98 %s %lf %lf \n", cAtomName, dnSLprot, dnSLdeut);
                
                fnCoord2ArraynSL(&iPosZ, dPosZ);				//calculate array pos
                //printf("Position z: %i dz: %lf", iPosZ, dPosZ);
                dMatrixnSL[iPosZ][0]+=dnSLprot;					//add nSL w/ exchange
                dMatrixnSL[iPosZ][1]+=dnSLdeut;					//add nSL w/o exchange
                i=iDimensionAtomNameTable;					    //end loop
                
            }
            else if ((strcmp(cAtomNameTable[i][2],"99") == 0) && (COMPLEX == 1))						//Non-exchangeable Hydrogen but complex of deuterated and non-deuterated protein
            {
                if (strcmp(cChainIdentifier, PROTCHAIN) == 0){
                    dnSLprot=dAtomPropertyTable[iDimensionAtomPropertyTable-1][1];
                }
                else {
                    dnSLprot=dAtomPropertyTable[iDimensionAtomPropertyTable-2][1];
                }
                dnSLdeut=dnSLprot;
                //printf("nSLFill98 %s %lf %lf \n", cAtomName, dnSLprot, dnSLdeut);
                
                fnCoord2ArraynSL(&iPosZ, dPosZ);				//calculate array pos
                //printf("Position z: %i dz: %lf", iPosZ, dPosZ);
                dMatrixnSL[iPosZ][0]+=dnSLprot;					//add nSL w/ exchange
                dMatrixnSL[iPosZ][1]+=dnSLdeut;					//add nSL w/o exchange
                i=iDimensionAtomNameTable;					    //end loop
                
            }
            else
            {
                iFoundInAtomPropertyTable=0;				     //different atom than hydrogen.
                i2=atoi(cAtomNameTable[i][2]);
                for (i3=0; i3<iDimensionAtomPropertyTable; i3++) //find element in property table
                    if (i2==dAtomPropertyTable[i3][0]) 
                    {
                        iFoundInAtomPropertyTable=1;
                        dnSLprot=dAtomPropertyTable[i3][1];			//look up nSL
                        dnSLdeut=dnSLprot;
                        //printf("nSLFill %s %lf %lf \n", cAtomName, dnSLprot, dnSLdeut);
                        
                        fnCoord2ArraynSL(&iPosZ, dPosZ);				//calculate array pos
                        dMatrixnSL[iPosZ][0]+=dnSLprot;					//add nSL w/ exchange
                        dMatrixnSL[iPosZ][1]+=dnSLdeut;					//add nSL w/o exchange
                        i=iDimensionAtomNameTable;					//end loop
                        
                    }
                if (iFoundInAtomPropertyTable==0) {printf("Unknown atom indentifier %i \n", i2);}
            }
        }
    }
    if (iFoundInAtomNameTable==0) 
    {
        printf("Unknown Atom %s %s \n", cResidue, cAtomName);
    }
    
    
}


//-------------------------------------------------------------------------------------------------------------------------


void fnFillSolvent(int iX, int iY, int iZ)
{
    
    int iNotYetFilled, iXl, iYl, iZl;
    
    //fill start cell with solvent
    iRecursiveCallCounter=0;				//initialize, for limiting recursive calls
    iSolventAtomsFilled=0;				//initialize, for progress monitoring
    iLastFraction=0;						//initialize, for progress monitoring
    
    fnFillProbeVolume(iX,iY,iZ);			//call recursive filling procedure
    
    iNotYetFilled=1;						//initialize
    
    while (iNotYetFilled==1)				//repeat as long as there are saved starting points
    {
        iNotYetFilled=0;					
        for (iXl=0; iXl<MATRIXSIZE; iXl++)			//search whole matrix
            for (iYl=0; iYl<MATRIXSIZE; iYl++)
                for (iZl=0; iZl<MATRIXSIZE; iZl++)
                {
                    if (iMatrix[iXl][iYl][iZl]==4)			//saved starting point?
                    {
                        iMatrix[iXl][iYl][iZl]=2;				//remove starting point label;
                        fnFillProbeVolume(iXl,iYl,iZl);		//call recursive procedure with this starting point 
                        iNotYetFilled=1;						//there might be other startingn points
                    }
                }
        fnAnalyzeMatrix();								//print intermediate results		  
    }
}

//-------------------------------------------------------------------------------------------------------------------------


void fnLoadFile(char* strFilename, int structure)								// load file and occupy matrix
{
    int iNumber, i, c, iNumberOfStrucInFile, iAtomCounter;
    char cAtomType[6], cChainString[10], cNumber[10];
    char cString[80], cAtomName[10], cResidue[10];
    double dDouble;
    double dPosX, dPosY, dPosZ;
    
    FILE *fp;									
    fp=fopen(strFilename,"r");
    
    iNumberOfStrucInFile=0;
    iAtomCounter=0;
    
    if(fp==NULL) {
        printf("Error: can't open file.\n");}
    
    
    while(!feof(fp)) {
        strcpy(cString,"");;
        fscanf(fp, "%3s", cString);
        if (strcmp(cString,"END") == 0) {
            iNumberOfStrucInFile+=1;
            printf("Increased iNumberOfStrucInFile: %i \n",iNumberOfStrucInFile);
            printf("\n");
            do
            {
                c = fgetc(fp);
            } while((c!='\n') && (c!=EOF));
        }
        else {
            fseek(fp, strlen(cString)*(-1), SEEK_CUR);                                       //rewind the 3 characters and evaluate
            fscanf(fp, "%6c", cString);
            for (i=0; i<6; i++) {
                if (cString[i]==' ') {
                    cString[i]='\0';
                }
            }
            cString[6]='\0';
            
            if (strcmp (cString,"ATOM") == 0)
            {
                
                fscanf(fp,"%6c", cNumber); //reading information
                cNumber[6]='\0';
                fscanf(fp,"%5c", cAtomName); //reading information
                cAtomName[5]='\0';
                while (cAtomName[0]==' ') {
                    for (i=0; i<6; i++) {
                        cAtomName[i]=cAtomName[i+1];
                    }
                }
                for (i=0; i<5; i++) {
                    if (cAtomName[i]==' ') {
                        cAtomName[i]='\0';
                    }
                }
                
                fscanf(fp,"%4s", cResidue); //reading information
                fscanf(fp,"%2s", cChainString); //reading information
                iNumber=atoi(cNumber);
                
                if ((cChainString[0]=='A')||(cChainString[0]=='B')||(cChainString[0]=='C')||(cChainString[0]=='D')||(cChainString[0]=='E')||(cChainString[0]=='F')||(cChainString[0]=='G')||(cChainString[0]=='H')||(cChainString[0]=='I')||(cChainString[0]=='J')||(cChainString[0]=='K')||(cChainString[0]=='L')||(cChainString[0]=='M')||(cChainString[0]=='N')||(cChainString[0]=='O')||(cChainString[0]=='P')||(cChainString[0]=='Q')||(cChainString[0]=='R')||(cChainString[0]=='S')||(cChainString[0]=='T')||(cChainString[0]=='U')||(cChainString[0]=='V')||(cChainString[0]=='W')||(cChainString[0]=='X')||(cChainString[0]=='Y')||(cChainString[0]=='Z')) {                                    //chain identifier present
                    fscanf(fp,"%i %lf %lf %lf %lf %lf %s",
                           &iNumber, &dPosX, &dPosY, &dPosZ, &dDouble, &dDouble, cAtomType);
                    iAtomCounter++;
                    //printf("%s %i \n",cAtomType, iAtomCounter);
                }
                else {                                         //chain identifier not present
                    fscanf(fp,"%lf %lf %lf %lf %lf %s",
                           &dPosX, &dPosY, &dPosZ, &dDouble, &dDouble, cAtomType);
                }
                
                //if (strlen(cAtomType)>2) {fseek(fp, (strlen(cAtomType)+1)*(-1), SEEK_CUR);}  //rewind file if Atom Type not present
                
                //get cAtomType from cAtomName if possible
                if (cAtomName[0]=='H') {strcpy(cAtomType,"H");}
                if (cAtomName[0]=='C') {strcpy(cAtomType,"C");}
                if (cAtomName[0]=='N') {strcpy(cAtomType,"N");}
                if (cAtomName[0]=='O') {strcpy(cAtomType,"O");}
                if (cAtomName[0]=='S') {strcpy(cAtomType,"S");}
                
                
                if (strcmp (cResidue,"HSP") == 0) {strcpy(cResidue,"HIS");}    //account for CHARMM residues
                if (strcmp (cResidue,"HSD") == 0) {strcpy(cResidue,"HIS");}
                if (strcmp (cResidue,"HSE") == 0) {
                    strcpy(cResidue,"HIS");
                }
                
                if (iNumberOfStrucInFile==structure) {
                    fnFillVanDerWaals(dPosX, dPosY, dPosZ, cAtomName, cResidue);  //filling volume
                    fnFillnSlProfile(dPosZ, cResidue, cAtomName, cChainString);                 //filling nSL profile
                    //printf("done.");
                }
                
            }
            else if (strcmp (cString,"HETATM") == 0)
            {
                fscanf(fp,"%5i %4s %3s %s", &iNumber, cAtomName, cResidue, cChainString); //reading information
                
                if ((cChainString[0]=='A')||(cChainString[0]=='B')||(cChainString[0]=='C')||(cChainString[0]=='D')||(cChainString[0]=='E')||(cChainString[0]=='F')||(cChainString[0]=='G')||(cChainString[0]=='H')||(cChainString[0]=='I')||(cChainString[0]=='J')||(cChainString[0]=='K')||(cChainString[0]=='L')||(cChainString[0]=='M')||(cChainString[0]=='N')||(cChainString[0]=='O')||(cChainString[0]=='P')||(cChainString[0]=='Q')||(cChainString[0]=='R')||(cChainString[0]=='S')||(cChainString[0]=='T')||(cChainString[0]=='U')||(cChainString[0]=='V')||(cChainString[0]=='W')||(cChainString[0]=='X')||(cChainString[0]=='Y')||(cChainString[0]=='Z')) {                                    //chain identifier present
                    fscanf(fp,"%i %lf %lf %lf %lf %lf %s",
                           &iNumber, &dPosX, &dPosY, &dPosZ, &dDouble, &dDouble, cAtomType);
                }
                else {                                         //chain identifier not present
                    fscanf(fp,"%lf %lf %lf %lf %lf %s",
                           &dPosX, &dPosY, &dPosZ, &dDouble, &dDouble, cAtomType);
                }
                
                //fscanf(fp,"%5i %4s %lf %lf %lf",
                //&iNumber, cAtomName, &dPosX, &dPosY, &dPosZ);              //reading information, more simple structure for HETATM
                
                strcpy(cResidue,"NON");
                
                if (iNumberOfStrucInFile==structure) {
                    fnFillnSlProfile(dPosZ, cResidue, cAtomType, cChainString);  //filling nSL profile
                    fnFillVanDerWaals(dPosX, dPosY, dPosZ, cAtomType, cResidue);  //filling volume
                }
            }
            do
            {
                c = fgetc(fp);
            } while((c!='\n') && (c!=EOF));
        }
    }
    
    fclose(fp);
    
};

//-------------------------------------------------------------------------------------------------------------------------


int fnDetermineLimits(char* strFilename)  //determining coordinate extrema and checking if array size is sufficient
{
    int iNumber,iResidueNumber,iFirstRead,i, iResidueNumberCounter, iResidueNumberOld, iNumberOfStrucInFile, c;
    char cAtomName[5], cAtomType[5];
    char cString[80], cResidue[5], cChainString[10];
    double d, dDouble;
    double dPosX, dPosY, dPosZ;
    int iAminoAcidCounter[20];
    double dAminoAcidFraction[20];
    
    for (i=0; i<20; i++) {iAminoAcidCounter[i]=0;}
    
    iResidueNumberCounter=0; //no residue read yet
    iResidueNumberOld=-1;
    iFirstRead=0;
    iNumberOfStrucInFile=0;
    
    FILE *fp;
    fp=fopen(strFilename,"r");
    
    if(fp==NULL) {
        printf("Error: can't open file.\n");}
    
    
    while(!feof(fp)) {
        strcpy(cString,"");;
        fscanf(fp, "%3s", cString);
        if (strcmp(cString,"END") == 0) {
            iNumberOfStrucInFile+=1;
            printf("Increased iNumberOfStrucInFile: %i \n",iNumberOfStrucInFile);
            printf("\n");
            do
            {
                c = fgetc(fp);
            } while((c!='\n') && (c!=EOF));
        }
        else {
            fseek(fp, strlen(cString)*(-1), SEEK_CUR);                                       //rewind the 3 characters and evaluate
            fscanf(fp, "%6s", cString);
            if (strcmp (cString,"ATOM") == 0)
            {
                
                fscanf(fp,"%5i %4s %3s %s", &iNumber, cAtomName, cResidue, cChainString);    //reading information
                
                if ((cChainString[0]=='A')||(cChainString[0]=='B')||(cChainString[0]=='C')||(cChainString[0]=='D')||(cChainString[0]=='E')||(cChainString[0]=='F')||(cChainString[0]=='G')||(cChainString[0]=='H')||(cChainString[0]=='I')||(cChainString[0]=='J')||(cChainString[0]=='K')||(cChainString[0]=='L')||(cChainString[0]=='M')||(cChainString[0]=='N')||(cChainString[0]=='O')||(cChainString[0]=='P')||(cChainString[0]=='Q')||(cChainString[0]=='R')||(cChainString[0]=='S')||(cChainString[0]=='T')||(cChainString[0]=='U')||(cChainString[0]=='V')||(cChainString[0]=='W')||(cChainString[0]=='X')||(cChainString[0]=='Y')||(cChainString[0]=='Z')) {                                    //chain identifier present
                    fscanf(fp,"%i %lf %lf %lf %lf %lf %s",
                           &iResidueNumber, &dPosX, &dPosY, &dPosZ, &dDouble, &dDouble, cAtomType);
                    //printf("ResidueNumber    : %lf \n", (double)iResidueNumber);

                }
                else {                                                                       //chain identifier not present
                    iResidueNumber=i;
                    fscanf(fp,"%lf %lf %lf %lf %lf %s",
                           &dPosX, &dPosY, &dPosZ, &dDouble, &dDouble, cAtomType);
                    
                }
                
                //if (strlen(cString)>2) {fseek(fp, strlen(cString)*(-1), SEEK_CUR);}  //rewind file if Atom Type not present
                
                if (strcmp (cResidue,"HSP") == 0) {strcpy(cResidue,"HIS");}    //account for CHARM residues
                if (strcmp (cResidue,"HSD") == 0) {strcpy(cResidue,"HIS");}
                if (strcmp (cResidue,"HSE") == 0) {
                    strcpy(cResidue,"HIS");
                }
                
                //first read coordinate sets initial limits
                if (iFirstRead==0) {dPosXMin=dPosX; dPosXMax=dPosX; dPosYMin=dPosY; dPosYMax=dPosY; dPosZMin=dPosZ; dPosZMax=dPosZ; iFirstRead=1;}
                
                if (dPosX<dPosXMin) dPosXMin=dPosX;		//Do new coordinates extend limits?
                if (dPosY<dPosYMin) dPosYMin=dPosY;
                if (dPosZ<dPosZMin) dPosZMin=dPosZ;
                if (dPosX>dPosXMax) dPosXMax=dPosX;
                if (dPosY>dPosYMax) dPosYMax=dPosY;
                if (dPosZ>dPosZMax) dPosZMax=dPosZ;
                
                if (iResidueNumber!=iResidueNumberOld)			//make statistic of residues for theoretical Volume and nSLD
                {
                    iResidueNumberCounter++;
                    iResidueNumberOld=iResidueNumber;
                    for (i=0; i<21; i++) {
                        if (i==21) {
                            printf("Cannot identify residue %s \n",cResidue);
                        }
                        if (strcmp (cResidue,cNameAminoAcid[i]) == 0) {
                            iAminoAcidCounter[i]++;
                            break;
                        }
                    }
                    //printf("ResidueNumber: %lf ResidueCounter: %lf Residue: %s \n", (double)iResidueNumber, (double)iResidueNumberCounter, cResidue);

                }
                
            }
            else if (strcmp (cString,"HETATM") == 0)
            {
                
                fscanf(fp,"%5i %4s %3s %s", &iNumber, cAtomName, cResidue, cChainString); //reading information
                
                if ((cChainString[0]=='A')||(cChainString[0]=='B')||(cChainString[0]=='C')||(cChainString[0]=='D')||(cChainString[0]=='E')||(cChainString[0]=='F')||(cChainString[0]=='G')||(cChainString[0]=='H')||(cChainString[0]=='I')||(cChainString[0]=='J')||(cChainString[0]=='K')||(cChainString[0]=='L')||(cChainString[0]=='M')||(cChainString[0]=='N')||(cChainString[0]=='O')||(cChainString[0]=='P')||(cChainString[0]=='Q')||(cChainString[0]=='R')||(cChainString[0]=='S')||(cChainString[0]=='T')||(cChainString[0]=='U')||(cChainString[0]=='V')||(cChainString[0]=='W')||(cChainString[0]=='X')||(cChainString[0]=='Y')||(cChainString[0]=='Z')) {                                    //chain identifier present
                    fscanf(fp,"%i %lf %lf %lf %lf %lf %s",
                           &iNumber, &dPosX, &dPosY, &dPosZ, &dDouble, &dDouble, cAtomType);
                }
                else {                                         //chain identifier not present
                    fscanf(fp,"%lf %lf %lf %lf %lf %s",
                           &dPosX, &dPosY, &dPosZ, &dDouble, &dDouble, cAtomType);
                }
                
                //fscanf(fp,"%5i %4s %lf %lf %lf",
                //&iNumber, cAtomName, &dPosX, &dPosY, &dPosZ);              //reading information, more simple structure for HETATM
                
                if (iFirstRead==0) {dPosXMin=dPosX; dPosXMax=dPosX; dPosYMin=dPosY; dPosYMax=dPosY; dPosZMin=dPosZ; dPosZMax=dPosZ; iFirstRead=1;}
                
                if (dPosX<dPosXMin) dPosXMin=dPosX;		//Do new coordinates extend limits?
                if (dPosY<dPosYMin) dPosYMin=dPosY;
                if (dPosZ<dPosZMin) dPosZMin=dPosZ;
                if (dPosX>dPosXMax) dPosXMax=dPosX;
                if (dPosY>dPosYMax) dPosYMax=dPosY;
                if (dPosZ>dPosZMax) dPosZMax=dPosZ;
            }
        }
    }
    
    fclose(fp);

    printf("Number of Structures in File: %i \n",iNumberOfStrucInFile);
    printf("Total Number of Residues per structure: %i \n",iResidueNumberCounter/iNumberOfStrucInFile);

    printf("Residue Summary \n \n");
    printf("RES \t # \t fraction \t nSLD \n");
    for (i=0; i<20; i++) 
    { 
        d=(double)iAminoAcidCounter[i];
        dDouble=(double)iResidueNumberCounter;
        dAminoAcidFraction[i]=d/dDouble;
        printf("%s \t %i \t %lf \t %e \n", cNameAminoAcid[i], iAminoAcidCounter[i]/iNumberOfStrucInFile, dAminoAcidFraction[i], dnSLDAminoAcid[i]);
    }
    
    printf("\n \n");
    
    d=0;
    for (i=0; i<20; i++) {d+=iAminoAcidCounter[i]*dnSLAminoAcid[i];} 
    printf("Total nSL of the protein expected from sequence [Å]: %e \n \n", d / (double)iNumberOfStrucInFile);
    
    d=0;
    for (i=0; i<20; i++) {d+=dAminoAcidFraction[i]*dnSLDAminoAcid[i];}
    printf("Average nSLD of the protein expected from sequence [Å^-2]: %e \n", d);
    
    d=0;
    for (i=0; i<20; i++) {d+=iAminoAcidCounter[i]*dVolumeAminoAcid[i];}
    printf("Volume of the protein expected from sequence [Å^3]: %e \n \n", d / (double)iNumberOfStrucInFile);
    
    printf("Size limits  \n xmin, xmax: %lf %lf \n ymin, ymax: %lf %lf \n zmin, zmax:  %lf %lf \n\n", 
           dPosXMin, dPosXMax, dPosYMin, dPosYMax, dPosZMin, dPosZMax);
    
    if (((dPosXMax-dPosXMin+2*iProbeFrame)>MATRIXSIZE*UNITSIZE) || 
        ((dPosYMax-dPosYMin+2*iProbeFrame)>MATRIXSIZE*UNITSIZE) || 
        ((dPosZMax-dPosZMin+2*iProbeFrame)>MATRIXSIZE*UNITSIZE))
    {
        printf("Matrix size too small. Change parameter MATRIXSIZE");
        return 1;
    }
    
    return 0;
}

//-----------------------------------------------------------------------------------------------------------------------
int strcpprot(char *cAtomName, char *cTableEntry)    //manually build string comparison for comparison of atom name with table
//returns 0 when entry found
{
    int i,j,k;
    char cCompareString[]="a";
    char cNumStr1[]="00000";
    char cNumStr2[]="00000";
    
    j=0;k=0;
    
    if (strcmp(cAtomName,cTableEntry)==0) {return 0;}     //identical identifiers -> positive return
    
    if (strlen(cAtomName)<=strlen(cTableEntry))			 //check if AtomName can be found in TableEntry
    {
        for (i=0; i<strlen(cAtomName); i++)					//process character by character
        {
            cCompareString[0]=cAtomName[i];
            if (!(strpbrk(cCompareString,cTableEntry)))		//character must be in table entry -> this sorts out wrong atoms
            {return 1;}
            
            if ((cAtomName[i]<58) && (cAtomName[i]>47))		//a number?
            {
                cNumStr1[j]=cAtomName[i];
                j++;
            }
            
            if ((cTableEntry[i]<58) && (cTableEntry[i]>47))		//a number?
            {
                cNumStr2[k]=cTableEntry[i];
                k++;
            }
        }
        //if there are two numbers in the cAtomName and cTableEntry they must 
        //be swapped difference between IUPAC and PDB nomenclature
        if (j!=k) {return 1;}                                         //number of numbers must be identical in both strings
        
        if (j==2) {
            i=cNumStr1[0]; cNumStr1[0]=cNumStr1[1]; cNumStr1[1]=i;
        }
		
        if (strcmp(cNumStr1,cNumStr2)==0) {
            //printf("%s \t %s \t %s \t %s \t %i \n", cAtomName, cTableEntry, cNumStr1, cNumStr2, j);
            return 0;
        }
		
        return 1;
    }
    
    return 1;
}
//-------------------------------------------------------------------------------------------------------------------------

void fnStoreResults(int structure)
{
    int iX, iY, iZ, iZnSL, iProteinCounter, iSolventCounter, iVoidCounter;
    double frac1, frac2;
    
    iZnSL=0;
    for(iZ=iProbeFrame; iZ<(MATRIXSIZE-iProbeFrame); iZ++)
    {
        iProteinCounter=0; iSolventCounter=0; iVoidCounter=0;
        for(iY=0; iY<MATRIXSIZE; iY++)
            for(iX=0; iX<MATRIXSIZE; iX++)
            {
                switch(iMatrix[iX][iY][iZ])
                {
                    case 0: iVoidCounter++; break;
                    case 1: iProteinCounter++; break;
                    case 2: iSolventCounter++; break;
                }
            }
        frac1=1/((double)(structure)+1);   //weighing factor for new results
        frac2=1-frac1;                     //weighing factor for old results
        dMatrixStore[iZnSL][0]=dMatrixStore[iZnSL][0]*frac2+dMatrixnSL[iZnSL][0]*frac1;
        dMatrixStore[iZnSL][1]=dMatrixStore[iZnSL][1]*frac2+dMatrixnSL[iZnSL][1]*frac1;
        dMatrixStore[iZnSL][2]=dMatrixStore[iZnSL][2]*frac2+(double)(iProteinCounter+iVoidCounter)*UNITSIZE*UNITSIZE*frac1;
        iZnSL+=1;
    }
    
}

//-------------------------------------------------------------------------------------------------------------------------

void fnWriteVolumeSliced(char* strFilename)	//write out results
{
    int iZ, iZnSL;
    
    FILE *fp;
    fp=fopen(strFilename,"w");
    
    if(fp==NULL) {printf("Error: can't write file.\n");}
    else {
    
        printf("Volumes written to file in Volume per delta z. This equals an area (of each slab). \n");
        printf("The solvent excluded volume in the output file contains the void volume and the protein volume. \n");
        
        fprintf(fp, "z   protnSL   deutnSL   area\n \n");
        iZnSL=0;
        for(iZ=iProbeFrame; iZ<(MATRIXSIZE-iProbeFrame); iZ++)		
        {	
            fprintf(fp, "%lf %lf %lf %lf \n", (double)iZ*UNITSIZE, dMatrixStore[iZnSL][0], dMatrixStore[iZnSL][1], dMatrixStore[iZnSL][2]);
            iZnSL+=1;
        }        
        
        fclose(fp);
        
    }
    
    
}

//-------------------------------------------------------------------------------------------------------------------------


int main (int argc, const char * argv[]) {
    
    int iX, iY, iZ, tilt, orientation, structure;
    char strFilename[300], buf[5];
    
    
    iProbeRadius=floor(PROBERADIUS/UNITSIZE+0.5);                   //integer probe radius used for distance comparisons on the grid 
    iProbeFrame=ceil(PROBERADIUS/UNITSIZE)+1;                       //the number of bins the grid is extended in each of the six directions
                                                                    //in order to have a complete solvation shell around the pdb structure.
                                                                    //this extension (frame) must have at least the size of half a solvent 
                                                                    //molecule 

    if (BATCHMODE==0){
        
        if (fnDetermineLimits(filename)==1) {return 0;}                         //check file and choosen array size
        
        printf("Loading file and occupying matrix ...\n\n");
        
        for (structure=0; structure<NUMBEROFSTRUCTURES; structure++) {
            
            for(iX=0; iX<MATRIXSIZE; iX++)                                         //initialize arrays
                for(iY=0; iY<MATRIXSIZE; iY++)
                    for(iZ=0; iZ<MATRIXSIZE; iZ++) {
                        iMatrix[iX][iY][iZ]=0;}
            
            for(iZ=0; iZ<MATRIXSIZE; iZ++) {
                dMatrixnSL[iZ][0]=0;
                dMatrixnSL[iZ][1]=0;
            }
            
            fnLoadFile(filename,structure);                                            //load certain structure to file and occupy array
            fnAnalyzeMatrix();                                                         //print intermediate results
            
            printf("Filling solvent volume for structure %i...\n\n", structure);
            fnFillSolvent(0,0,0);									                   //fill solvent volume
            fnAnalyzeMatrix();									                       //print intermediate results
            fnStoreResults(structure);
        }
        
        printf("Writing Output ...\n\n");						       //write output file
        fnWriteVolumeSliced(FilenameWrite);
    }
    else {
        
        for (tilt=TILTSTART; tilt<TILTEND; tilt+=TILTSTEP) {
            for (orientation=ORIESTART; orientation<ORIEEND; orientation+=ORIESTEP) {
                
                
                strcpy(strFilename,"");
                strcat(strFilename,strPath);
                strcat(strFilename,strFileNameRoot);
                strcat(strFilename,"_tilt");
                sprintf(buf,"%i",tilt);
                strcat(strFilename,buf);
                strcat(strFilename,"_orie");
                sprintf(buf,"%i",orientation);
                strcat(strFilename,buf);
                strcat(strFilename,".pdb");
                
                if (fnDetermineLimits(strFilename)==1) {return 0;}                         //check file and choosen array size
                
                for (structure=0; structure<NUMBEROFSTRUCTURES; structure++) {
                    
                    printf("Computing file %s...\n",strFilename);
                    printf("structure %i ...\n\n", structure);
                    
                    for(iX=0; iX<MATRIXSIZE; iX++)                                         //initialize arrays
                        for(iY=0; iY<MATRIXSIZE; iY++)
                            for(iZ=0; iZ<MATRIXSIZE; iZ++) {
                                iMatrix[iX][iY][iZ]=0;}
                    
                    for(iZ=0; iZ<MATRIXSIZE; iZ++) {
                        dMatrixnSL[iZ][0]=0;
                        dMatrixnSL[iZ][1]=0;
                    }                    
                    
                    fnLoadFile(strFilename,structure);                                         //load certain structure to file and occupy array
                    fnAnalyzeMatrix();                                                         //print intermediate results
                    
                    fnFillSolvent(0,0,0);									                   //fill solvent volume
                    fnAnalyzeMatrix();									                       //print intermediate results
                    fnStoreResults(structure);
                }                

                strcpy(strFilename,"");
                strcat(strFilename,strPath);
                strcat(strFilename,strFileNameRoot);
                strcat(strFilename,"_tilt");
                sprintf(buf,"%i",tilt);
                strcat(strFilename,buf);
                strcat(strFilename,"_orie");
                sprintf(buf,"%i",orientation);
                strcat(strFilename,buf);
                strcat(strFilename,".txt");

                printf("Writing Output ...\n\n");						       //write output file
                fnWriteVolumeSliced(strFilename);
            }
        }

        
    }
    
    
    return 0;
}
