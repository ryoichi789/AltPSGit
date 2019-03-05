#ifndef CBA_PDB_RECORD_H
#define CBA_PDB_RECORD_H

#include <string>
#include <vector>

namespace Cba {

typedef std::string PdbLine;

struct PdbRecord {
	PdbRecord();
	virtual ~PdbRecord();
	virtual void read(const PdbLine& line) = 0;

	void trunk_space(std::string& s);
};

// HEADER: Mandatory  single
struct PdbHeader : public PdbRecord {
	void read(const PdbLine& line);
	std::string classification;
	std::string depdate;
	std::string idcode;
};

// OBSLTE: Optional   single continued
struct PdbObslte : public PdbRecord {
	void read(const PdbLine& line);
	std::string repdate;
	std::string idcode;
	std::vector<std::string> ridcode;
};

// TITLE: Mandatory  single continued
struct PdbTitle : public PdbRecord {
	void read(const PdbLine& line);
	std::string title;
};

// CAVEAT: Optional   single continued
struct PdbCaveat : public PdbRecord {
	void read(const PdbLine& line);
	std::string idcode;
	std::string comment;
};

// COMPND: Mandatory  single continued
struct PdbCompnd : public PdbRecord {
	void read(const PdbLine& line);
	std::string compound;
};

// SOURCE: Mandatory  single continued
struct PdbSource : public PdbRecord {
	void read(const PdbLine& line);
	std::string source;
};

// KEYWDS: Mandatory  single continued
struct PdbKeywds : public PdbRecord {
	void read(const PdbLine& line);
	std::vector<std::string> keywords;
};

// EXPDTA: Mandatory  single continued
struct PdbExpdta : public PdbRecord {
	void read(const PdbLine& line);
	std::string technique;
};

// AUTHOR: Mandatory  single continued
struct PdbAuthor : public PdbRecord {
	void read(const PdbLine& line);
	std::vector<std::string> authors;
};

// REVDAT: Mandatory  multiple
struct PdbRevdat : public PdbRecord {
	void read(const PdbLine& line);
	int modnum;
	std::string moddate;
	std::string modid;
	int modtype;
	std::vector<std::string> record;
};

// SPRSDE: Optional   single continued
struct PdbSprsde : public PdbRecord {
	void read(const PdbLine& line);
	std::string sprsdedate;
	std::string idcode;
	std::vector<std::string> sidcode;
};

// fields used in JRNL and REMARK1 records
struct PdbReference {
	std::string auth;
	std::string titl;
	std::string edit;
	std::string ref;
	std::string publ;
	std::string refn;
};

// JRNL: Optional   other (single continued)
struct PdbJrnl : public PdbRecord {
	void read(const PdbLine& line);
	PdbReference reference;
};

// REMARK 1: Optional   other
struct PdbRemark1 : public PdbRecord {
	void read(const PdbLine& line);
	std::vector<PdbReference> references;
};


// REMARK 2: Mandatory  other
struct PdbRemark2 : public PdbRecord {
	void read(const PdbLine& line);
	float resolution;
	std::string comment;
};

// REMARK 3: Mandatory  other
struct PdbRemark3 : public PdbRecord {
	void read(const PdbLine& line);
	std::string info;
};

// REMARK N (N>=4): Optional   other
struct PdbRemarkN : public PdbRecord {
	void read(const PdbLine& line);
	int remarknum;
	std::string comment;
};

// DBREF: Optional   multiple
struct PdbDbref : public PdbRecord {
	void read(const PdbLine& line);
	std::string idcode;
	char chainid;
	int seqbegin;
	char insertbegin;
	int seqend;
	char insertend;
	std::string database;
	std::string dbaccession;
	std::string dbidcode;
	int dbseqbegin;
	char dbinsbeg;
	int dbseqend;
	char dbinsend;
};

// SEQADV: Optional   multiple
struct PdbSeqadv : public PdbRecord {
	void read(const PdbLine& line);
	std::string idcode;
	std::string resname;
	char chainid;
	int seqnum;
	char icode;
	std::string database;
	std::string dbidcode;
	std::string dbres;
	int dbseq;
	std::string conflict;
};

// SEQRES: Optional   multiple
struct PdbSeqres : public PdbRecord {
	void read(const PdbLine& line);
	char chainid;
	int numres;
	std::vector<std::string> resname;
};

// MODRES: Optional   multiple
struct PdbModres : public PdbRecord {
	void read(const PdbLine& line);
	std::string idcode;
	std::string resname;
	char chainid;
	int seqnum;
	char icode;
	std::string stdres;
	std::string comment;
};

// HET: Optional   multiple
struct PdbHet : public PdbRecord {
	void read(const PdbLine& line);
	std::string hetid;
	char chainid;
	int seqnum;
	char icode;
	int numhetatoms;
	std::string text;
};

// HETNAM: Optional   multiple continued
struct PdbHetnam : public PdbRecord {
	void read(const PdbLine& line);
	std::string hetid;
	std::string text;
};

// HETSYN: Optional   multiple
struct PdbHetsyn : public PdbRecord {
	void read(const PdbLine& line);
	std::string hetid;
	std::string hetsynonyms;
};

// FORMUL: Optional   multiple continued
struct PdbFormul : public PdbRecord {
	void read(const PdbLine& line);
	int compnum;
	std::string hetid;
	char asterisk;
	std::string text;
};

// HELIX: Optional   multiple
struct PdbHelix : public PdbRecord {
	void read(const PdbLine& line);
	enum Type {
		R_ALPHA = 1,
		R_OMEGA = 2,
		R_PI = 3,
		R_GAMMA = 4,
		R_310 = 5,
		L_ALPHA = 6,
		L_OMEGA = 7,
		L_GAMMA = 8,
		_27_RIBBON = 9,
		POLYPROLINE = 10
	};
	std::string technique;
	int sernum;
	std::string helixid;
	std::string initresname;
	char initchainid;
	int initseqnum;
	char initicode;
	std::string endresname;
	char endchainid;
	int endseqnum;
	char endicode;
	int helixclass;
	std::string comment;
	int length;
};

// SHEET: Optional   multiple
struct PdbSheet : public PdbRecord {
	void read(const PdbLine& line);
	int strand;
	std::string initresname;
	char initchainid;
	int initseqnum;
	char initicode;
	std::string endresname;
	char endchainid;
	int endseqnum;
	char endicode;
	int sense;
	std::string curatom;
	std::string curresname;
	char curchainid;
	int curresseq;
	char curicode;
	std::string prevatom;
	std::string prevresname;
	char prevchainid;
	int prevresseq;
	char previcode;
};

// TURN: Optional   multiple
struct PdbTurn : public PdbRecord {
	void read(const PdbLine& line);
	std::string turnid;
	std::string initresname;
	char initchainid;
	int initseqnum;
	char initicode;
	std::string endresname;
	char endchainid;
	int endseqnum;
	char endicode;
	std::string comment;
};

// SSBOND: Optional   multiple
struct PdbSsbond : public PdbRecord {
	void read(const PdbLine& line);
	char chainid1;
	int seqnum1;
	char icode1;
	char chainid2;
	int seqnum2;
	char icode2;
	std::string sym1;
	std::string sym2;
};

// LINK: Optional   multiple
struct PdbLink : public PdbRecord {
	void read(const PdbLine& line);
	std::string name1;
	char altloc1;
	std::string resname1;
	char chainid1;
	int resseq1;
	char icode1;
	std::string name2;
	char altloc2;
	std::string resname2;
	char chainid2;
	int resseq2;
	char icode2;
	std::string sym1;
	std::string sym2;
};

// HYDBND: Optional   multiple
struct PdbHydbnd : public PdbRecord {
	void read(const PdbLine& line);
	std::string name1;
	char altloc1;
	std::string resname1;
	char chainid1;
	int resseq1;
	char icode1;
	std::string nameh;
	char altloch;
	char chainh;
	int resseqh;
	char icodeh;
	std::string name2;
	char altloc2;
	std::string resname2;
	char chainid2;
	int resseq2;
	char icode2;
	std::string sym1;
	std::string sym2;
};

// SLTBRG: Optional   multiple
struct PdbSltbrg : public PdbRecord {
	void read(const PdbLine& line);
	std::string atom1;
	char altloc1;
	std::string resname1;
	char chainid1;
	int resseq1;
	char icode1;
	std::string atom2;
	char altloc2;
	std::string resname2;
	char chainid2;
	int resseq2;
	char icode2;
	std::string sym1;
	std::string sym2;
};

// CISPEP: Optional   multiple
struct PdbCispep : public PdbRecord {
	void read(const PdbLine& line);
	std::string pep1;
	char chainid1;
	int seqnum1;
	char icode1;
	std::string pep2;
	char chainid2;
	int seqnum2;
	char icode2;
	int modnum;
	float measure;
};

// SITE: Optional   multiple
struct PdbSite : public PdbRecord {
	void read(const PdbLine& line);
	int seqnum;
	std::string siteid;
	int numres;
	std::string resname1;
	char chainid1;
	int seq1;
	char icode1;
	std::string resname2;
	char chainid2;
	int seq2;
	char icode2;
	std::string resname3;
	char chainid3;
	int seq3;
	char icode3;
	std::string resname4;
	char chainid4;
	int seq4;
	char icode4;
};

// CRYST1: Mandatory  single
struct PdbCryst1 : public PdbRecord {
	void read(const PdbLine& line);
	float a;
	float b;
	float c;
	float alpha;
	float beta;
	float gamma;
	std::string sgroup;
	int z;
};

// ORIGX1 ORIGX2 ORIGX3: Mandatory  single
struct PdbOrigx : public PdbRecord {
	void read(const PdbLine& line);
	float o[3][3];
	float t[3];
	PdbOrigx();
	PdbOrigx(const PdbOrigx& rec);
	PdbOrigx& operator=(const PdbOrigx& rec);
};

// SCALE1 SCALE2 SCALE3: Mandatory  single
struct PdbScale : public PdbRecord {
	void read(const PdbLine& line);
	float s[3][3];
	float u[3];
	PdbScale();
	PdbScale(const PdbScale& rec);
	PdbScale& operator=(const PdbScale& rec);
};

// MTRIX1 MTRIX2 MTRIX3: Optional   multiple
struct PdbMtrix : public PdbRecord {
	void read(const PdbLine& line);
	int serial;
	float m[3][3];
	float v[3];
	int igiven;
	PdbMtrix();
	PdbMtrix(const PdbMtrix& rec);
	PdbMtrix& operator=(const PdbMtrix& rec);
};

// TVECT: Optional   multiple
struct PdbTvect : public PdbRecord {
	void read(const PdbLine& line);
	int serial;
	float t1;
	float t2;
	float t3;
	std::string text;
};

// MODEL: Optional   grouping
struct PdbModel : public PdbRecord {
	void read(const PdbLine& line);
	int serial;
};

struct PdbEndmdl : public PdbRecord {
	void read(const PdbLine& line);
};

// ATOM, HETATM: Optional   multiple
struct PdbAtom : public PdbRecord {
	PdbAtom() : hetatm(false) {}
	void read(const PdbLine& line);
	int serial;
	std::string name;
	char altloc;
	std::string resname;
	char chainid;
	int resseq;
	char icode;
	float x;
	float y;
	float z;
	float occupancy;
	float tempfactor;
	std::string segid;
	std::string element;
	std::string charge;
	bool hetatm; // true if HETATM
};

// SIGATM: Optional   multiple
struct PdbSigatm : public PdbRecord {
	void read(const PdbLine& line);
	int serial;
	std::string name;
	char altloc;
	std::string resname;
	char chainid;
	int resseq;
	char icode;
	float sigx;
	float sigy;
	float sigz;
	float sigocc;
	float sigtemp;
	std::string segid;
	std::string element;
	std::string charge;
};

// ANISOU: Optional   multiple
struct PdbAnisou : public PdbRecord {
	void read(const PdbLine& line);
	int serial;
	std::string name;
	char altloc;
	std::string resname;
	char chainid;
	int resseq;
	char icode;
	int u00;
	int u11;
	int u22;
	int u01;
	int u02;
	int u12;
	std::string segid;
	std::string element;
	std::string charge;
};

// SIGUIJ: Optional   multiple
struct PdbSiguij : public PdbRecord {
	void read(const PdbLine& line);
	int serial;
	std::string name;
	char altloc;
	std::string resname;
	char chainid;
	int resseq;
	char icode;
	int sig11;
	int sig22;
	int sig33;
	int sig12;
	int sig13;
	int sig23;
	std::string segid;
	std::string element;
	std::string charge;
};

// TER                     Optional   grouping
struct PdbTer : public PdbRecord {
	void read(const PdbLine& line);
	int serial;
	std::string resname;
	char chainid;
	int resseq;
	char icode;
};

// CONECT: Optional   multiple
struct PdbConect : public PdbRecord {
	void read(const PdbLine& line);
	int serial;
	int serial_bonded1;
	int serial_bonded2;
	int serial_bonded3;
	int serial_bonded4;
	int serial_hbonded1;
	int serial_hbonded2;
	int serial_saltbridged1;
	int serial_hbonded3;
	int serial_hbonded4;
	int serial_saltbridged2;
};

// MASTER: Mandatory  single
struct PdbMaster : public PdbRecord {
	void read(const PdbLine& line);
	int numremark;
	int numhet;
	int numhelix;
	int numsheet;
	int numturn;
	int numsite;
	int numxform;
	int numcoord;
	int numter;
	int numconect;
	int numseq;
};

// END: Mandatory  single
struct PdbEnd : public PdbRecord {
	void read(const PdbLine& line);
};

// unknown record
struct PdbUnknownRecord : public PdbRecord {
	void read(const PdbLine& line);
};

}
#endif
