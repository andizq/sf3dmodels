/*
 *  py_utils.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#ifndef PY_UTILS_H
#define PY_UTILS_H

#include <Python.h>

#include "inpars.h"

#define PY_MAX_LEN_PAR_NAME      30
#define PY_MAX_LEN_PAR_TYPE      20
#define PY_STR_LEN_0             80
#define UFUNC_BUFFER_SIZE       100

typedef struct {
  char name[PY_MAX_LEN_PAR_NAME+1], type[PY_MAX_LEN_PAR_TYPE+1];
  _Bool mandatory,isList;
} parTemplateType;

struct tempType{
  int intValue;
  double doubleValue;
  char strValue[PY_STR_LEN_0+1];
  _Bool boolValue,isNone;
};

extern PyObject *pMacros_global;

extern PyObject *pDensity,\
                *pTemperature,\
                *pAbundance,\
                *pMolNumDensity,\
                *pDoppler,\
                *pVelocity,\
                *pMagfield,\
                *pGasIIdust,\
                *pGridDensity;

extern PyObject *pModule_global;

extern _Bool userFuncsInitialized;

parTemplateType*	setTemplateDefaults(void);
void	myStrCpy(const char source[], char destination[], const int strlenDest);
int	getModuleFromName(char *moduleNameNoSuffix, PyObject **pModule);
int	getParTemplates(PyObject *pParsClassOrInstance, parTemplateType **parTemplates, int *nPars);
int	mallocInputParStrs(inputPars *par);
int	readParImg(PyObject *pPars, parTemplateType *parTemplates, const int nPars, parTemplateType *imgParTemplates, const int nImgPars, inputPars *par, image **img, int *nImages, void (*warning)(char *message));

int	setMacros(void);
void	unsetMacros(void);
void	getPythonFunc(PyObject *pModule, const char *funcName, PyObject **pFunc);
void	setUpUserPythonFuncs(PyObject *pModule);
void	setUpUserPythonFuncs_new(void);
void	decrefAllUserFuncs(void);
int	userFuncWrapper(PyObject *pFunc, const char *funcName\
  , double x, double y, double z, double *resultBuffer, int *numElemInUserFuncReturn);

void	pyFreeInputImgPars(image *inimg, int nImages);

#endif /* PY_UTILS_H */

