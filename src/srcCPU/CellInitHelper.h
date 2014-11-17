/*
 * CellInitHelper.h
 *
 *  Created on: Sep 22, 2013
 *      Author: wsun2
 */

#ifndef CELLINITHELPER_H_
#define CELLINITHELPER_H_

#include <vector>
#include "GeoVector.h"
#include <cmath>
#include "commonData.h"
#include "ConfigParser.h"
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <string>
#include "MeshGen.h"

using namespace std;

/**
 * Parameter the controls the simualtion.
 */
struct SimulationGlobalParameter {
public:
	std::string animationNameBase;
	double totalSimuTime;
	double dt;
	int totalTimeSteps;
	int totalFrameCount;
	int aniAuxVar;
	AnimationCriteria aniCri;
	std::string dataOutput;
	std::string imgOutput;
	std::string dataFolder;
	std::string dataName;
	void initFromConfig();
};

/**
 * Handles cell initialization.
 */
class CellInitHelper {
	SimulationType simuType;
	vector<CVector> internalBdryPts;

	CVector getPointGivenAngle(double currentAngle, double r,
			CVector centerPos);
	void generateRandomAngles(vector<double> &randomAngles,
			int initProfileNodeSize);

	void generateCellInitNodeInfo_v2(vector<CVector> &initPos);
	void generateECMInitNodeInfo(vector<CVector> &initECMNodePoss,
			int initNodeCountPerECM);
	void generateECMCenters(vector<CVector> &ECMCenters,
			vector<CVector> &CellCenters, vector<CVector> &bdryNodes);

	bool anyECMCenterTooClose(vector<CVector> &ecmCenters, CVector position);
	bool anyCellCenterTooClose(vector<CVector> &cellCenters, CVector position);
	bool anyBoundaryNodeTooClose(vector<CVector> &bdryNodes, CVector position);

	/**
	 * generate a random number between min and max.
	 */
	double getRandomNum(double min, double max);

	/**
	 * generates an array that could qualify for initial position of nodes in a cell.
	 */
	vector<CVector> generateInitCellNodes();

	/**
	 * Attempt to generate an array that represents relative position of nodes in a cell.
	 */
	vector<CVector> attemptGeenerateInitCellNodes();

	/**
	 * Determine if an array, representing relative position of nodes in a cell,
	 * is qualified for initialization purpose.
	 */
	bool isPositionQualify(vector<CVector> &poss);

	/**
	 * Initialize internal boundary points given input file.
	 */
	void initInternalBdry();

	/**
	 * Initialize raw input given an array of cell center positions.
	 */
	void initializeRawInput(RawDataInput &rawInput,
			std::vector<CVector> &cellCenterPoss);

	/**
	 * generate initial node positions given cartilage raw data.
	 */
	void transformRawCartData(CartilageRawData &cartRawData, CartPara &cartPara,
			std::vector<CVector> &initNodePos);

	/**
	 * Given a cell center position, decide if the center position is MX type or not.
	 */
	bool isMXType(CVector position);

	vector<CVector> rotate2D(vector<CVector> &initECMNodePoss, double angle);

	/**
	 * Used for generate RawDataInput for actual simulation purpose.
	 * The raw data generated does not only contain cell center positions but also
	 * epithelium node positions and boundary node positions.
	 */
	RawDataInput generateRawInputWithProfile(
			std::vector<CVector> &cellCenterPoss, bool isInnerBdryIncluded =
					true);

	/**
	 * generate simulation initialization data _v2 given raw data.
	 */
	SimulationInitData_V2 initInputsV3(RawDataInput &rawData);

	/**
	 * Used for generate RawDataInput for stabilization purpose.
	 * In order to place cells in a random domain, we need an initialization step to stabilize
	 * cell center positions which are generated by triangular mesh. It is usually difficult for
	 * triangular meshing to generate semi-evenly displaced cell centers.
	 */
	RawDataInput generateRawInput_stab();

	/**
	 * Used for generate RawDataInput for single cell model test purpose.
	 */
	RawDataInput generateRawInput_singleCell();

	/**
	 * Used for generate RawDataInput for actual simulation purpose.
	 */
	RawDataInput generateRawInput_simu(std::vector<CVector> &cellCenterPoss);

	/**
	 * generate simulation initialization data given raw data.
	 */
	SimulationInitData initInputsV2(RawDataInput &rawData);

public:
	/**
	 * Default constructor is used.
	 */
	CellInitHelper();
	/**
	 * Default de-constructor is used.
	 */
	virtual ~CellInitHelper();

	/**
	 * generate simulation initialization data _v2.
	 */
	SimulationInitData_V2 initStabInput();

	/**
	 * generate simulation initialization data _v2 given raw data.
	 */
	SimulationInitData_V2 initSimuInput(std::vector<CVector> &cellCenterPoss);

	/**
	 * generate simulation initialization data _v2 for single cell testing.
	 */
	SimulationInitData_V2 initSingleCellTest();
};

#endif /* CELLINITHELPER_H_ */
