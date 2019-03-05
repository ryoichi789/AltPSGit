#include <iostream>
#include <string>
#include <cstdlib>
#include <cctype>
#include "pdb_reader.h"
#include "pdb_entry.h"
#include "pdb_record.h"

using namespace Cba;
using std::string;

// read:
bool PdbReader::read(PdbEntry& ent, std::istream& from)
{
	ent.clear();

	PdbLine line;
	int l;
	while (std::getline(from, line)) {
		if ((l = line.length()) < 80) line.append(80-l,' ');
		add_record(ent, line);
	}
	ent.make_entry();

	if (ent.nrecords()) return true;
	return false;
}

// get_record:
void PdbReader::add_record(PdbEntry& ent, const std::string& line)
{
	// add a new record
	// unless line is continued from previous line
	if (line.at(0) == ' ') {
		return;
	} else if (line.compare(0,6,"HEADER") == 0) {
		ent.m_records.push_back(ent.m_header = new PdbHeader);
	} else if (line.compare(0,6,"OBSLTE") == 0) {
		if (line.at(9) == ' ')
			ent.m_records.push_back(ent.m_obslte = new PdbObslte);
	} else if (line.compare(0,6,"TITLE ") == 0) {
		if (line.at(9) == ' ')
			ent.m_records.push_back(ent.m_title = new PdbTitle);
	} else if (line.compare(0,6,"CAVEAT") == 0) {
		if (line.at(9) == ' ')
			ent.m_records.push_back(ent.m_caveat = new PdbCaveat);
	} else if (line.compare(0,6,"COMPND") == 0) {
		if (line.at(9) == ' ')
			ent.m_records.push_back(ent.m_compnd = new PdbCompnd);
	} else if (line.compare(0,6,"SOURCE") == 0) {
		if (line.at(9) == ' ')
			ent.m_records.push_back(ent.m_source = new PdbSource);
	} else if (line.compare(0,6,"KEYWDS") == 0) {
		if (line.at(9) == ' ')
			ent.m_records.push_back(ent.m_keywds = new PdbKeywds);
	} else if (line.compare(0,6,"EXPDTA") == 0) {
		if (line.at(9) == ' ')
			ent.m_records.push_back(ent.m_expdta = new PdbExpdta);
	} else if (line.compare(0,6,"AUTHOR") == 0) {
		if (line.at(9) == ' ')
			ent.m_records.push_back(ent.m_author = new PdbAuthor);
	} else if (line.compare(0,6,"REVDAT") == 0) {
		if (line.at(11) == ' ') {
			ent.m_revdats.push_back(new PdbRevdat);
			ent.m_records.push_back(ent.m_revdats.back());
		}
	} else if (line.compare(0,6,"SPRSDE") == 0) {
		if (line.at(9) == ' ')
			ent.m_records.push_back(ent.m_sprsde = new PdbSprsde);
	} else if (line.compare(0,6,"JRNL  ") == 0) {
		if (line.at(17) == ' ')
			ent.m_records.push_back(ent.m_jrnl = new PdbJrnl);
	} else if (line.compare(0,10,"REMARK   1") == 0) {
		if (ent.m_records.empty() ||
			!dynamic_cast<PdbRemark1*>(ent.m_records.back()))
			ent.m_records.push_back(ent.m_remark1 = new PdbRemark1);
	} else if (line.compare(0,10,"REMARK   2") == 0) {
		if (ent.m_records.empty() ||
			!dynamic_cast<PdbRemark2*>(ent.m_records.back()))
			ent.m_records.push_back(ent.m_remark2 = new PdbRemark2);
	} else if (line.compare(0,10,"REMARK   3") == 0) {
		if (ent.m_records.empty() ||
			!dynamic_cast<PdbRemark3*>(ent.m_records.back()))
			ent.m_records.push_back(ent.m_remark3 = new PdbRemark3);
	} else if (line.compare(0,6,"REMARK") == 0) {
		ent.m_remarkns.push_back(new PdbRemarkN);
		ent.m_records.push_back(ent.m_remarkns.back());
		//PdbRemarkN* r = dynamic_cast<PdbRemarkN*>(ent.m_records.back());
		//if (!r || (r->remarknum != std::atoi(line.substr(7,3).c_str()))) {
		//	ent.m_remarkns.push_back(new PdbRemarkN);
		//	ent.m_records.push_back(ent.m_remarkns.back());
		//}
	} else if (line.compare(0,6,"DBREF ") == 0) {
		ent.m_dbrefs.push_back(new PdbDbref);
		ent.m_records.push_back(ent.m_dbrefs.back());
	} else if (line.compare(0,6,"SEQADV") == 0) {
		ent.m_seqadvs.push_back(new PdbSeqadv);
		ent.m_records.push_back(ent.m_seqadvs.back());
	} else if (line.compare(0,6,"SEQRES") == 0) {
		PdbSeqres* r;
		if (ent.m_records.empty() ||
			!(r = dynamic_cast<PdbSeqres*>(ent.m_records.back())) ||
			r->chainid != line.at(11)) {
			ent.m_seqreses.push_back(new PdbSeqres);
			ent.m_records.push_back(ent.m_seqreses.back());
		}
	} else if (line.compare(0,6,"MODRES") == 0) {
		ent.m_modreses.push_back(new PdbModres);
		ent.m_records.push_back(ent.m_modreses.back());
	} else if (line.compare(0,6,"HET   ") == 0) {
		ent.m_hets.push_back(new PdbHet);
		ent.m_records.push_back(ent.m_hets.back());
	} else if (line.compare(0,6,"HETNAM") == 0) {
		if (line.substr(8,2).compare("  ") == 0) {
			ent.m_hetnams.push_back(new PdbHetnam);
			ent.m_records.push_back(ent.m_hetnams.back());
		}
	} else if (line.compare(0,6,"HETSYN") == 0) {
		if (line.at(9) == ' ') {
			ent.m_hetsyns.push_back(new PdbHetsyn);
			ent.m_records.push_back(ent.m_hetsyns.back());
		}
	} else if (line.compare(0,6,"FORMUL") == 0) {
		if (line.at(17) == ' ') {
			ent.m_formuls.push_back(new PdbFormul);
			ent.m_records.push_back(ent.m_formuls.back());
		}
	} else if (line.compare(0,6,"HELIX") == 0) {
		ent.m_helices.push_back(new PdbHelix);
		ent.m_records.push_back(ent.m_helices.back());
	} else if (line.compare(0,6,"SHEET") == 0) {
		ent.m_sheets.push_back(new PdbSheet);
		ent.m_records.push_back(ent.m_sheets.back());
	} else if (line.compare(0,6,"TURN") == 0) {
		ent.m_turns.push_back(new PdbTurn);
		ent.m_records.push_back(ent.m_turns.back());
	} else if (line.compare(0,6,"SSBOND") == 0) {
		ent.m_ssbonds.push_back(new PdbSsbond);
		ent.m_records.push_back(ent.m_ssbonds.back());
	} else if (line.compare(0,6,"LINK  ") == 0) {
		ent.m_links.push_back(new PdbLink);
		ent.m_records.push_back(ent.m_links.back());
	} else if (line.compare(0,6,"HYDBND") == 0) {
		ent.m_hydbnds.push_back(new PdbHydbnd);
		ent.m_records.push_back(ent.m_hydbnds.back());
	} else if (line.compare(0,6,"SLTBRG") == 0) {
		ent.m_sltbrgs.push_back(new PdbSltbrg);
		ent.m_records.push_back(ent.m_sltbrgs.back());
	} else if (line.compare(0,6,"CISPEP") == 0) {
		ent.m_cispeps.push_back(new PdbCispep);
		ent.m_records.push_back(ent.m_cispeps.back());
	} else if (line.compare(0,6,"SITE  ") == 0) {
		ent.m_sites.push_back(new PdbSite);
		ent.m_records.push_back(ent.m_sites.back());
	} else if (line.compare(0,6,"CRYST1") == 0) {
		ent.m_records.push_back(ent.m_cryst1 = new PdbCryst1);
	} else if (line.compare(0,5,"ORIGX") == 0) {
		if (line.at(5) == '1')
			ent.m_records.push_back(ent.m_origx = new PdbOrigx);
	} else if (line.compare(0,5,"SCALE") == 0) {
		if (line.at(5) == '1')
			ent.m_records.push_back(ent.m_scale = new PdbScale);
	} else if (line.compare(0,5,"MTRIX") == 0) {
		if (line.at(5) == '1') {
			ent.m_mtrices.push_back(new PdbMtrix);
			ent.m_records.push_back(ent.m_mtrices.back());
		}
	} else if (line.compare(0,6,"TVECT ") == 0) {
		ent.m_tvects.push_back(new PdbTvect);
		ent.m_records.push_back(ent.m_tvects.back());
	} else if (line.compare(0,6,"MODEL ") == 0) {
		//if (line.compare(10, 4, "   0") == 0) return; // can happen, e.g., in 1CTL
		ent.m_models.push_back(new PdbModel);
		ent.m_records.push_back(ent.m_models.back());
	} else if (line.compare(0,6,"ENDMDL") == 0) {
		ent.m_endmdls.push_back(new PdbEndmdl);
		ent.m_records.push_back(ent.m_endmdls.back());
	} else if (line.compare(0,6,"ATOM  ") == 0) {
		if (line.at(16) == 'B' || line.at(16) == '2') return;
		ent.m_atoms.push_back(new PdbAtom);
		ent.m_records.push_back(ent.m_atoms.back());
	} else if (line.compare(0,6,"HETATM") == 0) {
		if (line.at(16) == 'B' || line.at(16) == '2') return;
		PdbAtom* r = new PdbAtom;
		r->hetatm = true;
		ent.m_atoms.push_back(r);
		ent.m_records.push_back(r);
	} else if (line.compare(0,6,"SIGATM") == 0) {
		ent.m_sigatms.push_back(new PdbSigatm);
		ent.m_records.push_back(ent.m_sigatms.back());
	} else if (line.compare(0,6,"ANISOU") == 0) {
		if (line.at(16) == 'B' || line.at(16) == '2') return;
		ent.m_anisous.push_back(new PdbAnisou);
		ent.m_records.push_back(ent.m_anisous.back());
	} else if (line.compare(0,6,"SIGUIJ") == 0) {
		ent.m_siguijs.push_back(new PdbSiguij);
		ent.m_records.push_back(ent.m_siguijs.back());
	} else if (line.compare(0,6,"TER   ") == 0) {
		ent.m_ters.push_back(new PdbTer);
		ent.m_records.push_back(ent.m_ters.back());
	} else if (line.compare(0,6,"CONECT") == 0) {
		ent.m_conects.push_back(new PdbConect);
		ent.m_records.push_back(ent.m_conects.back());
	} else if (line.compare(0,6,"MASTER") == 0) {
		ent.m_records.push_back(ent.m_master = new PdbMaster);
	} else if (line.compare(0,6,"END   ") == 0) {
		ent.m_records.push_back(ent.m_end = new PdbEnd);
	} else {
		ent.m_unknowns.push_back(new PdbUnknownRecord);
		ent.m_records.push_back(ent.m_unknowns.back());
	}

	ent.m_records.back()->read(line);
}
