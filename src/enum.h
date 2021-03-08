#ifndef ENUM_H
#define ENUM_H

#include <string>
#include <iostream>

using namespace std;

enum euvIdentifier {W171, W193, W195, W284, W304, WUNKNOWN};
typedef enum euvIdentifier euvID;

enum methodIdentifier {METH_NUMERIC, METH_SHC, METH_ELLIPTICAL, METH_UNKNOWN};
typedef enum methodIdentifier methodID;

enum modelIdentifier {MODEL_PFSS, MODEL_CSSS, MODEL_UNKNOWN};
typedef enum modelIdentifier modelID;

enum groupIdentifier {GROUP_STANFORD, GROUP_KIEL, GROUP_BALA, GROUP_UNKNOWN};
typedef enum groupIdentifier groupID;

enum spacecraftIdentifier {SC_ACE, SC_SOHO, SC_STEREOA, SC_STEREOB, SC_ULYSSES, SC_HELIOS, SC_VOYAGER1, SC_VOYAGER2, SC_NONE};
typedef enum spacecraftIdentifier spacecraftID;

enum instrumentIdentifier {	INST_NONE, INST_SOHO_LASCOC1, INST_SOHO_LASCOC2, INST_SOHO_LASCOC3,
							INST_SOHO_EIT171, INST_SOHO_EIT195, INST_SOHO_EIT284, INST_SOHO_EIT304,
							INST_ULYSSES_SWICS, INST_ULYSSES_VHMFGM, INST_ULYSSES_SWOOPS,
							INST_ACE_SWEPAM_full, INST_ACE_MAG4MIN, INST_ACE_SWESWI,
							INST_STEREO_PLASTIC, INST_STEREO_MAG1MIN};
typedef enum instrumentIdentifier instrumentID;

enum originIdentifier {ORIGIN_WSO, ORIGIN_NSOGONG, ORIGIN_NSOKPVT, ORIGIN_SOHOMDI, ORIGIN_SOHOMDI_DAILY, ORIGIN_VER_DAILY,
	ORIGIN_SDOHMI, ORIGIN_OWN, ORIGIN_UNKNOWN};
typedef enum originIdentifier originID;		// TODO origin of photospheric magnetograms, rename this struct

inline string getStringFromEUVid(euvID id)
{
	return   id == W171 ? "W171"	:
			(id == W193 ? "W193"	:
			(id == W195 ? "W195"	:
			(id == W284 ? "W284"	:
			(id == W304 ? "W304"	: "UnknownWL"))));
}

inline string getStringFromGroupID(groupID group)
{
	return 		group 	== GROUP_STANFORD 	? "Stanford"		:
			(	group 	== GROUP_KIEL		? "Kiel"			:
			(	group 	== GROUP_BALA		? "Bala"			: "UnknownGroup"));
}

inline string getStringFromModelID(modelID model)
{
	return		model 	== MODEL_PFSS		? "PFSS"			:
			(	model 	== MODEL_CSSS 		? "CSSS" 			: "UnknownModel");
}

inline string getStringFromMethodID(methodID method)
{
	return 		method 	== METH_NUMERIC		? "GridSph"			:
			(	method 	== METH_ELLIPTICAL	? "GridEll"			:
			(	method 	== METH_SHC			? "SHC"				: "UnknownMethod"));
}

inline string getStringFromSpacecraftID(spacecraftID id)
{
	return		id == SC_ACE				? "ACE"				:
			(	id == SC_HELIOS				? "HELIOS"			:
			(	id == SC_SOHO				? "SOHO"			:
			(	id == SC_STEREOA			? "STEREO-A"		:
			(	id == SC_STEREOB			? "STEREO-B"		:
			(	id == SC_ULYSSES			? "ULYSSES"			:
			(	id == SC_VOYAGER1			? "VOYAGER1"		:
			(	id == SC_VOYAGER1			? "VOYAGER2"		:
			(	id == SC_NONE				? "SC_NONE"			: "UnknownSpacecraft"))))))));
}

inline string getStringFromInstrumentID(instrumentID id)
{
	return 		id==INST_SOHO_LASCOC1		? "LASCO_C1"		:
			(	id==INST_SOHO_LASCOC2		? "LASCO_C2"		:
			(	id==INST_SOHO_LASCOC3		? "LASCO_C3"		:
			(	id==INST_SOHO_EIT171 		? "EIT_171" 		:
			(	id==INST_SOHO_EIT195 		? "EIT_195" 		:
			(	id==INST_SOHO_EIT284 		? "EIT_284" 		:
			(	id==INST_SOHO_EIT304 		? "EIT_304" 		:
			(	id==INST_ULYSSES_SWICS 		? "ULYSSES-SWICS" 	:
			(	id==INST_ULYSSES_SWOOPS		? "ULYSSES-SWOOPS"	:
			(	id==INST_ULYSSES_VHMFGM 	? "ULYSSES-VHMFGM" 	:
			(	id==INST_ACE_SWEPAM_full 	? "ACE-SWEPAM_full" :
			(	id==INST_ACE_SWESWI		 	? "ACE-SWESWI"		:
			(	id==INST_ACE_MAG4MIN		? "ACE-MAG_4min"	:
			(	id==INST_STEREO_MAG1MIN		? "STEREO-MAG_1min"	:
			(	id==INST_STEREO_PLASTIC		? "STEREO-PLASTIC"	: "UnknownInstrument"))))))))))))));
}

inline string getStringFromOriginID(originID id)
{
	return		id == ORIGIN_WSO			? "WSO"				:
			(	id == ORIGIN_NSOGONG		? "GONG"			:
			(	id == ORIGIN_NSOKPVT		? "KPVT"			:
			(	id == ORIGIN_SOHOMDI		? "MDI"				:
			(	id == ORIGIN_SOHOMDI_DAILY	? "MDIDAILY"		:
			(	id == ORIGIN_VER_DAILY		? "VERDAILY"		:
			(	id == ORIGIN_SDOHMI			? "HMI"				:
			(	id == ORIGIN_OWN			? "OWN"				: "UnkownOrigin")))))));
}



inline euvID getEUVidFromString(const string &str)
{
	return 	 str == "W171" ? W171 :
			(str == "W193" ? W193 :
			(str == "W195" ? W195 :
			(str == "W284" ? W284 :
			(str == "W304" ? W304 : WUNKNOWN))));
}

inline groupID getGroupIDfromString(const string &str)
{
	return 		str 	== "Stanford"		? GROUP_STANFORD	:
			(	str 	== "Kiel"			? GROUP_KIEL		:
			(	str 	== "Bala"			? GROUP_BALA		: GROUP_UNKNOWN));
}

inline modelID getModelIDfromString(const string &str)
{
	return		str 	== "PFSS"			? MODEL_PFSS		:
			(	str 	== "CSSS"			? MODEL_CSSS		: MODEL_UNKNOWN);
}

inline methodID getMethodIDfromString(const string &str)
{
	return 		str 	== "GridSph" 	|| str 	== "GridSpheric" 	? METH_NUMERIC		:
			(	str 	== "GridEll"	|| str 	== "GridElliptic"	? METH_ELLIPTICAL	:
			(	str 	== "SHC"			? METH_SHC			:	METH_UNKNOWN));
}

inline spacecraftID getSpacecrraftIDFromString(const string &sc)
{
	return		sc == "ACE"			? SC_ACE		:
			(	sc == "HELIOS"		? SC_HELIOS		:
			(	sc == "SOHO"		? SC_SOHO		:
			(	sc == "STEREO-A"	? SC_STEREOA	:
			(	sc == "STEREO-B"	? SC_STEREOB	:
			(	sc == "ULYSSES"		? SC_ULYSSES	:
			(	sc == "VOYAGER1"	? SC_VOYAGER1	:
			(	sc == "VOYAGER2"	? SC_VOYAGER2	:
			(	sc == "SC_NONE"		? SC_NONE		: SC_NONE))))))));
}

inline instrumentID getInstrumentIDfromString(const string &instr)
{
	return 		instr=="LASCO_C1"			? INST_SOHO_LASCOC1			:
			(	instr=="LASCO_C2"			? INST_SOHO_LASCOC2			:
			(	instr=="LASCO_C3"			? INST_SOHO_LASCOC3			:
			(	instr=="EIT_171" 			? INST_SOHO_EIT171 			:
			(	instr=="EIT_195" 			? INST_SOHO_EIT195 			:
			(	instr=="EIT_284" 			? INST_SOHO_EIT284 			:
			(	instr=="EIT_304" 			? INST_SOHO_EIT304 			:
			(	instr=="ULYSSES-SWICS"		? INST_ULYSSES_SWICS 		:
			(	instr=="ULYSSES-SWOOPS"		? INST_ULYSSES_SWOOPS		:
			(	instr=="ULYSSES-VHMFGM" 	? INST_ULYSSES_VHMFGM 		:
			(	instr=="ACE-SWEPAM_full"	? INST_ACE_SWEPAM_full 		:
			(	instr=="ACE-SWESWI"			? INST_ACE_SWESWI 			:
			(	instr=="ACE-MAG_4min"		? INST_ACE_MAG4MIN			:
			(	instr=="STEREO-MAG_1min"	? INST_STEREO_MAG1MIN		:
			(	instr=="STEREO-PLASTIC"		? INST_STEREO_PLASTIC		: INST_NONE))))))))))))));
}

inline originID getOriginIDfromString(const string &str)
{
	return		str	== "WSO"				? ORIGIN_WSO				:
			(	str == "GONG"				? ORIGIN_NSOGONG			:
			(	str == "KPVT"				? ORIGIN_NSOKPVT			:
			(	str == "MDI"				? ORIGIN_SOHOMDI			:
			(	str == "MDIDAILY"			? ORIGIN_SOHOMDI_DAILY		:
			(	str == "VERDAILY"			? ORIGIN_VER_DAILY			:
			(	str == "HMI"				? ORIGIN_SDOHMI				: ORIGIN_UNKNOWN))))));
}

#endif
