// MoleculeGasSim.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define Pointcount 3
#define BoundingX 1 //Größe in X Richtung
#define BoundingY 1 //GrößE IN y rICHTUNG

typedef struct {
	double X;
	double Y;
}Vector;

typedef struct {

	Vector Position;
	Vector Velocity;
	double dMass;
	double dRadius;

}Molecule;

Molecule CreateMolcule(double PosX, double PosY, double VelX, double VelY, double Mas, double Rad) { //Initializes Molecule

	Vector Pos;
	Pos.X = PosX;
	Pos.Y = PosY;

	Vector Vel;
	Vel.X = VelX;
	Vel.Y = VelY;

	Molecule res;
	res.Position = Pos;
	res.Velocity = Vel;
	res.dMass = Mas;
	res.dRadius = Rad;

	return res;
}



double SizeVector(const Vector * V1) { //Returns the Size of a Vector

	return sqrt(pow(V1->X,2) + pow(V1->Y,2));
};

double ScalarProduct(const Vector * V1, const Vector * V2) { //Calculates Scalar Product of two Vectors

	return V1->X * V2->X + V1->Y * V2->Y;
}

Vector SubtractVector(const Vector * V1, const Vector * V2) { //Subtracts V2 from V1
	Vector result;
	result.X = V1->X - V2->X;
	result.Y = V1->Y - V2->Y;

	return result;

}

Vector NormVector(const Vector * V1) {
	Vector vRes;
	double dSize = SizeVector(V1);
	vRes.X = V1->X / dSize;
	vRes.Y = V1->Y / dSize;
	return vRes;
}

Vector SwapToNewBasis(const Vector * swapVector, const Vector * Basis1, const Vector * Basis2) { //Swap a Vector from Standard Basis to new one
	Vector vRes;
	vRes.X = ScalarProduct(swapVector, Basis1);
	vRes.Y = ScalarProduct(swapVector, Basis2);
	return vRes;
}

Vector SwapBackBasis(const Vector * swapVector, const Vector * Basis1, const Vector * Basis2) { //Swaps Vector Back to standard Basis
	Vector vRes;
	vRes.X = swapVector->X * Basis1->X + swapVector->Y *  Basis2->X;
	vRes.Y = swapVector->X * Basis1->Y + swapVector->Y *  Basis2->Y;
	return vRes;
}

Vector createVector(double X, double Y) { //constructor for a vector
	Vector res;
	res.X = X;
	res.Y = Y;
	return res;
}

double CalcDT(double d, double p1x, double p1y, double v1x, double v1y, double p2x, double p2y, double v2x, double v2y) { //calculate T from Input

	return (p1x*v1x - p1x * v2x - p2x * v1x + p2x * v2x + p1y * v1y - p1y * v2y - p2y * v1y + p2y * v2y + sqrt(pow(d, 2)* pow(v1x, 2) - 2 * pow(d, 2)*v1x*v2x + pow(d, 2)* pow(v2x, 2) + pow(d, 2)* pow(v1y, 2) - 2 * pow(d, 2)*v1y*v2y + pow(d, 2)* pow(v2y, 2) - pow(p1x, 2)* pow(v1y, 2) + 2 * pow(p1x, 2)*v1y*v2y - pow(p1x, 2)* pow(v2y, 2) + 2 * p1x*p2x*pow(v1y, 2) - 4 * p1x*p2x*v1y*v2y + 2 * p1x*p2x*pow(v2y, 2) + 2 * p1x*p1y*v1x*v1y - 2 * p1x*p1y*v1x*v2y - 2 * p1x*p1y*v2x*v1y + 2 * p1x*p1y*v2x*v2y - 2 * p1x*p2y*v1x*v1y + 2 * p1x*p2y*v1x*v2y + 2 * p1x*p2y*v2x*v1y - 2 * p1x*p2y*v2x*v2y - pow(p2x, 2)* pow(v1y, 2) + 2 * pow(p2x, 2)*v1y*v2y - pow(p2x, 2)* pow(v2y, 2) - 2 * p2x*p1y*v1x*v1y + 2 * p2x*p1y*v1x*v2y + 2 * p2x*p1y*v2x*v1y - 2 * p2x*p1y*v2x*v2y + 2 * p2x*p2y*v1x*v1y - 2 * p2x*p2y*v1x*v2y - 2 * p2x*p2y*v2x*v1y + 2 * p2x*p2y*v2x*v2y - pow(p1y, 2)* pow(v1x, 2) + 2 * pow(p1y, 2)*v1x*v2x - pow(p1y, 2)* pow(v2x, 2) + 2 * p1y*p2y*pow(v1x, 2) - 4 * p1y*p2y*v1x*v2x + 2 * p1y*p2y*pow(v2x, 2) - pow(p2y, 2)* pow(v1x, 2) + 2 * pow(p2y, 2)*v1x*v2x - pow(p2y, 2)* pow(v2x, 2))) / (pow(v1x, 2) - 2 * v1x*v2x + pow(v2x, 2) + pow(v1y, 2) - 2 * v1y*v2y + pow(v2y, 2));
}

void StepMolecule(Molecule * M1, double dt) { // Step the molecule time
	M1->Position.X = M1->Position.X + dt * M1->Velocity.X;
	M1->Position.Y = M1->Position.Y + dt * M1->Velocity.Y;
};

void MoleculeCollision(Molecule * M1, Molecule * M2) { //calculates collision betwean Molecules, with time reveresal this sould be exception free unles molecules have the exact same position and velocity

	double dBackstep; //Time we must go back to propperly calculate collisions

	dBackstep = CalcDT(M1->dRadius + M2->dRadius,M1->Position.X,M1->Position.Y,M1->Velocity.X,M1->Velocity.Y, M2->Position.X, M2->Position.Y, M2->Velocity.X, M2->Velocity.Y); //calc backwards timestep

	StepMolecule(M1, -1 * dBackstep); //"reverse" the time of the colliding molecules to the moment of collision
	StepMolecule(M2, -1 * dBackstep);

	//Calculete Colision, and !!New Velocity""
	
	Vector vNormal;	 //Vector connecting the two Masses
	Vector vOrthogonal;  //Vector orthogonal to connection

	 //Calcualte vNormal and Orthogonal
	{  
		//Create Vectors
		vNormal = SubtractVector(&M1->Position, &M2->Position);
		vOrthogonal.X = vNormal.Y;
		vOrthogonal.Y = -1 * vNormal.X;

		//Norm them
		vNormal = NormVector(&vNormal);
		vOrthogonal = NormVector(&vOrthogonal);
	}  //Calcualte vNormal and Orthogonal
	
	Vector vVelocity1NB, vVelocity2NB; //Current Velocity in New Basis
	Vector NewVel1NB, NewVel2NB; 

	 //Calcualte vVelocity1NB and 2NB
	{
		vVelocity1NB = SwapToNewBasis(&M1->Velocity, &vNormal, &vOrthogonal);
		vVelocity2NB = SwapToNewBasis(&M2->Velocity, &vNormal, &vOrthogonal);
	}

	//Calculate New Velocitys, Y component remains constant
	NewVel1NB.X = (M1->dMass * vVelocity1NB.X + M2->dMass *(2* vVelocity2NB.X - vVelocity1NB.X)) / (M1->dMass + M2->dMass);
	NewVel2NB.X = (M2->dMass * vVelocity2NB.X + M1->dMass *(2 * vVelocity1NB.X - vVelocity2NB.X)) / (M1->dMass + M2->dMass);
	
	NewVel1NB.Y = vVelocity1NB.Y;
	NewVel2NB.Y = vVelocity2NB.Y;


	M1->Velocity = SwapBackBasis(&NewVel1NB, &vNormal, &vOrthogonal);
	M2->Velocity = SwapBackBasis(&NewVel2NB, &vNormal, &vOrthogonal);

	StepMolecule(M1, dBackstep); //"walk" the time we reversed again
	StepMolecule(M2, dBackstep);

};

//Wall Collisions
void StiffCollision(Vector vNormal,Vector vWall , Molecule * M1) { //Collision with infinit Mass Object with normal Vector vNormal, vNormal must point into the volume and vWall is one point on the Wall !!!!vNormal musst be normed
	
	double dBackstep;

	if (1 != SizeVector(&vNormal)) { //Check if vNormal is Size 1 //TODO Maybe disable check for better performance? (Maybe check at compile time)
		printf("Normal Vector isnt normed!!!");
		exit(2); //Error vector wassnt normed
	}

	dBackstep = (vNormal.X*(M1->Position.X - vWall.X) - M1->dRadius + vNormal.Y*(M1->Position.Y - vWall.Y)) / (vNormal.X * M1->Velocity.X + vNormal.Y * M1->Velocity.Y);
	
	StepMolecule(M1, -1 * dBackstep); //"reverse" the time of the colliding molecules to the moment of collision


	Vector vOrthogonal;  //Vector orthogonal to connection

	 //Calcualte  vOrthogonal
	{
	
		vOrthogonal.X = vNormal.Y;
		vOrthogonal.Y = -1 * vNormal.X;

	}  //Calcualte vOrthogonal

	Vector vVelocity1NB; //Current Velocity in New Basis


	//Calcualte vVelocity1NB and 2NB
	{
		vVelocity1NB = SwapToNewBasis(&M1->Velocity, &vNormal, &vOrthogonal);
	}

	//Calculate New Velocitys, Y component remains constant
	
	vVelocity1NB.X = -1 * vVelocity1NB.X; //Invert Normal Speed

	M1->Velocity = SwapBackBasis(&vVelocity1NB, &vNormal, &vOrthogonal);
	

	StepMolecule(M1, dBackstep); //"walk" the time we reversed again

}

bool WallCollision(Molecule * M1) {
	if (M1->Position.X <= 0) { //Collision with left Wall
		StiffCollision(createVector(1, 0), createVector(0, 0), M1);
		return true;
	}
	else if (M1->Position.X >= BoundingX) { //Collision with right Wall
		StiffCollision(createVector(-1, 0), createVector(BoundingX, BoundingY), M1);
		return true;
	}
	else if (M1->Position.Y <= 0) { //Collision with lower Wall
		StiffCollision(createVector(0, 1), createVector(0,0), M1);
		return true;
	}
	else if (M1->Position.Y >= BoundingY) { //Collision with right Wall
		StiffCollision(createVector(0, -1), createVector(BoundingX, BoundingY), M1);
		return true;
	}
	return false; //Didnt Collide with a Wall
}

//*WallCollisions

Molecule Moc[Pointcount]; //hold data of all molecules that are part of the sim, is global for easy of access
int iTimestep; //Integer that holds the current timestep

void TimestepAll(double dT) { //Timesteps all molecules by dT

	for (int i = 0; i < Pointcount; ++i) {
		StepMolecule(Moc + i, dT);
	}

}

bool CheckCollision(Molecule * M1, Molecule * M2) { //Checks if two molecules have / are colliding
	Vector vNormal = SubtractVector(&M1->Position, &M2->Position);
	if (SizeVector(&vNormal) <= M1->dRadius + M2->dRadius) { //they are colliding
		return true;
	}
	return false; //they arent
}

void HandleCollisions(void) { //Handles collisions for all Molecules

	for (int i = 0; i < Pointcount; ++i) { //Run for all Molecules, with time revesal we musst allways check with all molecules
		bool collided;
		do {
			collided = false; //reset collision
			for (int j = 0; j < Pointcount; ++j) {
				if (j != i) { //Prevent selfe collision and division by 0
					if (CheckCollision(Moc + 1, Moc + j)) {
						MoleculeCollision(Moc + 1, Moc + j); //They have collided calc their behaviour and rerun loop
						collided = true; //They have collided, the loop must be rerun
						break;
					}
				}
				
			}
			collided = WallCollision(Moc + 1); //TODO Not optimal runns this check twice if at least one collision happens and then no collision happens, but seams acceptable
		} while (collided); //Rerun Loop if they have collided
	}

}

void WriteToCSV(void) {
	FILE * fp;
	char filename[50]; //String to hold the filename
	sprintf_s(filename, 50, "stepdata%d.csv", iTimestep); //Build Name of new file
	if (NULL != fopen_s(&fp, filename, "w")) {//open the file
		printf("FEHLER !!");
		exit(1); //Exit Code 1, couldnt open file
	}
	fprintf(fp, "X ; Y ; Vx ; Vy ;\n"); //File Header
	for (int i = 0; i < Pointcount; ++i) {
		fprintf(fp, "%f ; %f ; %f ; %f ;\n", Moc[i].Position.X, Moc[i].Position.Y, Moc[i].Velocity.X, Moc[i].Velocity.Y);
	}
	fclose(fp); //close File
	++iTimestep;
}

void FixedStepLoop(int Steps, double dT) {
	WriteToCSV(); //Save initial Data
	for (int i = 0; i < Steps; ++i) {
		TimestepAll(dT); //Advance Time
		HandleCollisions(); //Check and correct Collisions
		WriteToCSV(); //Save Data
	}
};

void initializeParticlesBasic(Vector vV, double dT, double dCp, double dMass, double dRadius, double dDensityDifference ) {
	//Initializes A box with a dDensityDiffrence chance to be located in the Left half with base velocity vV and Temp dependent velocity dT
	srand(time(NULL)); //Seed randome generator

	double dPx, dPy, dVx, dVy; //double to hold randome position

	double dET = sqrt(2 * dT*dCp);

	for (int i = 0; i < Pointcount; ++i) {
		
		dPx = ((double)rand()) / ((double)RAND_MAX);
		if (dPx <= dDensityDifference / (dDensityDifference + 1)) {
			dPx = 0.5 * dPx / (dDensityDifference / (dDensityDifference + 1));
		}else {
			dPx = 0.5f * (dPx- dDensityDifference / (dDensityDifference + 1)) / (1 / (dDensityDifference + 1)) + 0.5f;
		}
		dPx = dPx * BoundingX;
		dPy = ((double)rand()) / ((double)RAND_MAX) * BoundingY;

		dVx = (2 * ((double)rand()) / ((double)RAND_MAX) - 1) * dET;
		dVy = sqrt(pow(dET, 2) - pow(dVx, 2));
		if (((double)rand()) / ((double)RAND_MAX) > 0.5f) {
			dVy = -1 * dVy;
		}
		dVx = dVx + vV.X;
		dVy = dVy + vV.Y;

		Moc[i] = CreateMolcule(dPx, dPy, dVx, dVy, dMass, dRadius);

	}

}

int main()
{
	iTimestep = 0; //Set timestep to 0,
	initializeParticlesBasic(createVector(0,0),237,0.0001,1,0.001,2);
	FixedStepLoop(1, 0.001);

	return 0;
}

// Programm ausführen: STRG+F5 oder "Debuggen" > Menü "Ohne Debuggen starten"
// Programm debuggen: F5 oder "Debuggen" > Menü "Debuggen starten"

// Tipps für den Einstieg: 
//   1. Verwenden Sie das Projektmappen-Explorer-Fenster zum Hinzufügen/Verwalten von Dateien.
//   2. Verwenden Sie das Team Explorer-Fenster zum Herstellen einer Verbindung mit der Quellcodeverwaltung.
//   3. Verwenden Sie das Ausgabefenster, um die Buildausgabe und andere Nachrichten anzuzeigen.
//   4. Verwenden Sie das Fenster "Fehlerliste", um Fehler anzuzeigen.
//   5. Wechseln Sie zu "Projekt" > "Neues Element hinzufügen", um neue Codedateien zu erstellen, bzw. zu "Projekt" > "Vorhandenes Element hinzufügen", um dem Projekt vorhandene Codedateien hinzuzufügen.
//   6. Um dieses Projekt später erneut zu öffnen, wechseln Sie zu "Datei" > "Öffnen" > "Projekt", und wählen Sie die SLN-Datei aus.
