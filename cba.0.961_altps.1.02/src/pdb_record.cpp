#include <string>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include "pdb_record.h"

using namespace Cba;

//
PdbRecord::PdbRecord() {}

//
PdbRecord::~PdbRecord() {}

//
void PdbRecord::trunk_space(std::string& s)
{
	for (int i = 0; ; ) {
		i = s.find("  ", i);
		if (i >= s.length()) break;
		s.erase(i, 1);
	}
}

//
void PdbHeader::read(const PdbLine& line)
{
	classification = line.substr(10, 40);
	depdate = line.substr(50, 9);
	idcode = line.substr(62, 4);
}

//
void PdbTitle::read(const PdbLine& line)
{
	title.append(line, 10, 60);
	trunk_space(title);
}

//
void PdbCompnd::read(const PdbLine& line)
{
	compound.append(line, 10, 60);
	trunk_space(compound);
}

//
void PdbRemark2::read(const PdbLine& line)
{
	if (line.compare(11, 11, "RESOLUTION.") == 0 && std::isdigit(line.at(23)))
		resolution = std::atof(line.substr(22, 5).c_str());
	else
		resolution = 0.0;
}

//
void PdbSeqres::read(const PdbLine& line)
{
	chainid = line.at(11);
	numres = std::atoi(line.substr(13, 4).c_str());
	std::string a;
	for (int i = 19; i < 71; i += 4) {
		a = line.substr(i, 3);
		if (a.compare("   ") == 0) break;
		resname.push_back(a);
	}
}

//
void PdbHet::read(const PdbLine& line)
{
	hetid = line.substr(7, 3);
	chainid = line.at(12);
	seqnum = std::atoi(line.substr(13, 4).c_str());
	icode = line.at(17);
	numhetatoms = std::atoi(line.substr(20, 5).c_str());
	text = line.substr(30, 40);
	trunk_space(text);
	if (text.at(0) == ' ') text.erase(0, 1);
}

//
void PdbModel::read(const PdbLine& line)
{
	serial = std::atoi(line.substr(10, 4).c_str());
}

//
void PdbEndmdl::read(const PdbLine& line) {}

//
void PdbAtom::read(const PdbLine& line)
{
	serial = std::atoi(line.substr(6, 5).c_str());
	name = line.substr(12, 4);
	altloc = line.at(16);
	resname = line.substr(17, 3);
	chainid = line.at(21);
	resseq = std::atoi(line.substr(22, 4).c_str());
	icode = line.at(26);
	x = std::atof(line.substr(30, 8).c_str());
	y = std::atof(line.substr(38, 8).c_str());
	z = std::atof(line.substr(46, 8).c_str());
	occupancy = std::atof(line.substr(54, 6).c_str());
	tempfactor = std::atof(line.substr(60, 6).c_str());
	segid = line.substr(72, 4);
	element = line.substr(76, 2);
	charge = line.substr(78, 2);
	if (line.compare(0, 6, "HETATM") == 0) hetatm = true;
	else hetatm = false;
}

//
void PdbTer::read(const PdbLine& line) {}

//
void PdbEnd::read(const PdbLine& line) {}

//
void PdbObslte::read(const PdbLine& line) { /* to be implemented */ }
void PdbCaveat::read(const PdbLine& line) { /* to be implemented */ }
void PdbSource::read(const PdbLine& line) { /* to be implemented */ }
void PdbKeywds::read(const PdbLine& line) { /* to be implemented */ }
void PdbExpdta::read(const PdbLine& line) { /* to be implemented */ }
void PdbAuthor::read(const PdbLine& line) { /* to be implemented */ }
void PdbRevdat::read(const PdbLine& line) { /* to be implemented */ }
void PdbSprsde::read(const PdbLine& line) { /* to be implemented */ }
void PdbJrnl::read(const PdbLine& line) { /* to be implemented */ }
void PdbRemark1::read(const PdbLine& line) { /* to be implemented */ }
void PdbRemark3::read(const PdbLine& line) { /* to be implemented */ }
void PdbRemarkN::read(const PdbLine& line)
{
	//remarknum = std::atoi(line.substr(7,3).c_str());
	/* to be implemented */
}
void PdbDbref::read(const PdbLine& line) { /* to be implemented */ }
void PdbSeqadv::read(const PdbLine& line) { /* to be implemented */ }
void PdbModres::read(const PdbLine& line) { /* to be implemented */ }

//
void PdbHetnam::read(const PdbLine& line)
{
	hetid = line.substr(11, 3);
	text.append(line, 15, 55);
	trunk_space(text);
	if (text.at(0) == ' ') text.erase(0, 1);
}

//
void PdbHetsyn::read(const PdbLine& line)
{
	hetid = line.substr(11, 3);
	hetsynonyms.append(line, 15, 55);
	trunk_space(hetsynonyms);
	if (hetsynonyms.at(0) == ' ') hetsynonyms.erase(0, 1);
}

void PdbFormul::read(const PdbLine& line) { /* to be implemented */ }
void PdbHelix::read(const PdbLine& line) { /* to be implemented */ }
void PdbSheet::read(const PdbLine& line) { /* to be implemented */ }
void PdbTurn::read(const PdbLine& line) { /* to be implemented */ }
void PdbSsbond::read(const PdbLine& line) { /* to be implemented */ }
void PdbLink::read(const PdbLine& line) { /* to be implemented */ }
void PdbHydbnd::read(const PdbLine& line) { /* to be implemented */ }
void PdbSltbrg::read(const PdbLine& line) { /* to be implemented */ }
void PdbCispep::read(const PdbLine& line) { /* to be implemented */ }
void PdbSite::read(const PdbLine& line) { /* to be implemented */ }
void PdbCryst1::read(const PdbLine& line) { /* to be implemented */ }
void PdbOrigx::read(const PdbLine& line) { /* to be implemented */ }
void PdbScale::read(const PdbLine& line) { /* to be implemented */ }
void PdbMtrix::read(const PdbLine& line) { /* to be implemented */ }
void PdbTvect::read(const PdbLine& line) { /* to be implemented */ }
void PdbSigatm::read(const PdbLine& line) { /* to be implemented */ }
void PdbAnisou::read(const PdbLine& line) { /* to be implemented */ }
void PdbSiguij::read(const PdbLine& line) { /* to be implemented */ }
void PdbConect::read(const PdbLine& line) { /* to be implemented */ }
void PdbMaster::read(const PdbLine& line) { /* to be implemented */ }
void PdbUnknownRecord::read(const PdbLine& line) { /* do nothing */ }

//
PdbOrigx::PdbOrigx()
{
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) o[i][j] = 0;
	for (int i = 0; i < 3; i++) t[i] = 0;
}

//
PdbOrigx::PdbOrigx(const PdbOrigx& rec)
{
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) o[i][j] = rec.o[i][j];
	for (int i = 0; i < 3; i++) t[i] = rec.t[i];
}

//
PdbOrigx& PdbOrigx::operator=(const PdbOrigx& rec)
{
	if (this != &rec) {
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) o[i][j] = rec.o[i][j];
		for (int i = 0; i < 3; i++) t[i] = rec.t[i];
	}
	return *this;
}

//
PdbScale::PdbScale()
{
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) s[i][j] = 0;
	for (int i = 0; i < 3; i++) u[i] = 0;
}

//
PdbScale::PdbScale(const PdbScale& rec)
{
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) s[i][j] = rec.s[i][j];
	for (int i = 0; i < 3; i++) u[i] = rec.u[i];
}

//
PdbScale& PdbScale::operator=(const PdbScale& rec)
{
	if (this != &rec) {
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) s[i][j] = rec.s[i][j];
		for (int i = 0; i < 3; i++) u[i] = rec.u[i];
	}
	return *this;
}

//
PdbMtrix::PdbMtrix()
{
	serial = 0;
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) m[i][j] = 0;
	for (int i = 0; i < 3; i++) v[i] = 0;
	igiven = 0;
}

//
PdbMtrix::PdbMtrix(const PdbMtrix& rec)
{
	serial = rec.serial;
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) m[i][j] = rec.m[i][j];
	for (int i = 0; i < 3; i++) v[i] = rec.v[i];
	igiven = rec.igiven;
}

//
PdbMtrix& PdbMtrix::operator=(const PdbMtrix& rec)
{
	if (this != &rec) {
		serial = rec.serial;
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) m[i][j] = rec.m[i][j];
		for (int i = 0; i < 3; i++) v[i] = rec.v[i];
		igiven = rec.igiven;
	}
	return *this;
}
