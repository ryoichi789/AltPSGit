//  follow the NOTEs below
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
// NOTE: include necessary headers
#include "altpsm.h"
#include "sdf_entry.h"
#include "sdf_reader.h"
#include "pdb_entry.h"
#include "pdb_reader.h"
#include "molecule.h"
#include "protein.h"
#include <dirent.h>
#include <iomanip>


using namespace Cba;

AltPS the_app;

AltPS::AltPS() 
{
}

void AltPS::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tq: query PDB filename" << endl;
	cerr << "\tt: target PDB filename" << endl;
	cerr << "\to: output directory (current directory by default)" << endl;
	cerr << "\tsize: minimun size (atom number) for matching surface area [int]"
		<< " (30 by default)" << endl;
	cerr << "\tzscore: minimun z-score for matching surface area [float]"
		<< " (4.0 by default)" << endl;
}

void AltPS::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "q") == 0) {
		m_obj->m_infile1 = param.value;
	} else if (strcmp(param.key, "t") == 0) {
		m_obj->m_infile2 = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_dirname = param.value;
	} else if (strcmp(param.key, "size") == 0) {
		m_obj->m_cluster_natom = std::atoi(param.value);
	} else if (strcmp(param.key, "zscore") == 0) {
		m_obj->m_zscore = std::atof(param.value);
	} else {
		std::cerr << "unrecognized parameter "
			<< param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

void AltPS::set_param_default()
{
	m_prog = "altps";
	m_infile1 = "";
	m_infile2 = "";
	m_pairname = "";
	m_dirname = "./";
	m_score_smap = 0.80;
	m_score_threshold = 0.80;
	m_list_pair_maximum = 10;
	m_cluster_minimum = 3;
	m_cluster_natom = 30;
	m_zscore = 4.0;
}

void AltPS::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
}

int AltPS::run()
{
	std::cout << "query  pdbfile = " << m_infile1 << std::endl;
	std::cout << "target pdbfile = " << m_infile2 << std::endl;
	if(check_names() < 0) return -1;
	pairname(m_infile1,m_infile2);
	
	//compare file1 with file2
	std::ifstream q_ifs(m_infile1.c_str());
	PdbEntry q_ent;
	PdbReader::read(q_ent, q_ifs);
	Protein q_prot;
	for (int i = 0; i < q_ent.nmolecules(); i++) {
		if (q_ent.type_of_molecule(i) == PdbEntry::PROTEIN &&
			q_ent.create_protein(i, q_prot)) {
			m_mol1_id = q_prot.id();
			m_mol1_chainid = chain_id(q_prot);
			m_resnumbers_q = get_serial_number_of_residue(q_prot);
			std::vector<LocArea> q_lareas;
			std::vector<int> q_suratoms;
			m_dmat_q = calc_locarea(q_prot,q_lareas,q_suratoms);
			compare_with_pdb(q_prot,q_lareas,q_suratoms);
		}
	}
	return 0;
}

//check_names
int AltPS::check_names()
{
	int select = 1;
	if(m_infile1 == ""){
		select = -1;
		usage();
		std::cerr << "Please type query pdbfile name [-q]" << std::endl;
		return select;
	} 
	if(m_infile2 == ""){
		select == -1;
		usage();
		std::cerr << "Please type target pdbfile name [-t]" << std::endl;
	}
	std::ifstream fi1;
	fi1.open(m_infile1.c_str());
	if(!fi1){
		std::cerr << "Failed to open -q pdbfile" << std::endl;
	}
	std::ifstream fi2;
	fi2.open(m_infile2.c_str());
	if(!fi2){
		select = -1;
		std::cerr << "Failed to open -t pdbfile" << std::endl;
	}
	DIR *dp;
	if ( (dp=opendir(m_dirname.c_str())) == NULL ){
		select = -1;
		std::cerr << "Failed to open directory (-o directory)" << std::endl;
	}else{
		if(m_dirname.compare(m_dirname.size()-1,1,"/") != 0) m_dirname.append("/");
	}
	std::cout << "output directory = " << m_dirname << std::endl;
	return select;
}

// pairname
void AltPS::pairname(const std::string f1,const std::string f2)
{
	std::string fname1 = f1;
	std::string fname2 = f2;
	int num;
	if(fname1.find_last_of("/") != std::string::npos){
		num = fname1.find_last_of("/");fname1.erase(0,num+1);
	}
	if(fname1.find_last_of(".") != std::string::npos){
		num = fname1.find_last_of(".");fname1.erase(num);
	}
	if(fname2.find_last_of("/") != std::string::npos){
		num = fname2.find_last_of("/");fname2.erase(0,num+1);
	}
	if(fname2.find_last_of(".") != std::string::npos){
		num = fname2.find_last_of(".");fname2.erase(num);
	}
	m_pairname = fname1 + "_" + fname2;
}

// calc_locarea()
std::vector<std::vector<double> > AltPS::calc_locarea(Protein& prot, std::vector<LocArea>& lareas, std::vector<int>& suratoms)
{
	prot.assign_atom_types();
	prot.calc_atom_asas();
	std::vector<std::vector<double> > dmat = calc_dmat(prot);
	set_local_area(prot,lareas,suratoms,dmat);
	assign_property_larea(prot,lareas,suratoms);
	return dmat;
}

//calc_dmat
std::vector<std::vector<double> > AltPS::calc_dmat(Protein& prot)
{
	std::vector<std::vector<double> > dmat(prot.natoms());
	for(int i=0;i<prot.natoms();i++){ dmat[i].resize(prot.natoms(),0.0);}
	for (int i = 0; i < prot.natoms()-1; i++) {
	if(prot.atom(i)->asa > 0.0 ){
		for (int j = i+1; j < prot.natoms(); j++) {
		if(prot.atom(j)->asa > 0.0 ){
			double d = prot.atom(i)->pos.distance_from(prot.atom(j)->pos);
			dmat[i][j] = d;
			dmat[j][i] = d;
		}}
	}}
	return dmat;
}

//compare_with_pdb
void AltPS::compare_with_pdb(Protein& q_prot, std::vector<LocArea>& q_lareas, std::vector<int>& q_suratoms)
{
	std::ifstream t_ifs(m_infile2.c_str());
	PdbEntry t_ent;
	PdbReader::read(t_ent, t_ifs);
	Protein t_prot;
	for (int i = 0; i < t_ent.nmolecules(); i++) {
		if (t_ent.type_of_molecule(i) == PdbEntry::PROTEIN &&
			t_ent.create_protein(i, t_prot)) {
			m_mol2_id = t_prot.id();
			m_mol2_chainid = chain_id(t_prot);
			m_resnumbers_t = get_serial_number_of_residue(t_prot);
			std::vector<LocArea> t_lareas;
			std::vector<int> t_suratoms;
			m_dmat_t = calc_locarea(t_prot,t_lareas,t_suratoms);
			search_area(q_prot,q_lareas,q_suratoms,t_prot,t_lareas,t_suratoms);
			clustering(q_lareas,t_lareas);
			analyze_result(q_prot,q_lareas,t_prot,t_lareas);
		}
	}
}

//search_area
void AltPS::search_area(
	const Protein& q_prot, std::vector<LocArea>& q_lareas, std::vector<int>& q_suratoms,
	const Protein& t_prot, std::vector<LocArea>& t_lareas, std::vector<int>& t_suratoms)
{
	// get pair of local areas with TC3 >= m_score_smap : slist
	std::map<int, std::vector<int> > smap;
	for(int i = 0; i < q_suratoms.size(); i++){
		if(q_lareas[q_suratoms[i]].dense){
			std::vector<int> slist;
			std::vector<double> sslist;
			for(int j = 0; j < t_suratoms.size(); j++){
				if(t_lareas[t_suratoms[j]].dense && 
					patty_match( (q_prot.atom(q_suratoms[i]))->pc_class,
								 (t_prot.atom(t_suratoms[j]))->pc_class ))
				{
					double ss = q_lareas[q_suratoms[i]].tc3(t_lareas[t_suratoms[j]]);
					if(ss >= m_score_smap){
						slist.push_back(t_suratoms[j]);
						sslist.push_back(ss);
					}
				}
			}//j
			if(slist.size() > 0){
				if(slist.size() > 50) slist = get_vector_top(50,slist,sslist);
				smap.insert( std::map<int, std::vector<int> >::value_type(q_suratoms[i],slist));
			}
		}
	}
	// calculate TC1 score between atom i of local area q and atom j of t : m_qt_sco
	m_qt_sco.clear();
	for(int i = 0; i < q_suratoms.size(); i++){
		for(int j = 0; j < t_suratoms.size(); j++){
			std::vector<int> qt_pair;
			qt_pair.push_back(q_suratoms[i]);
			qt_pair.push_back(t_suratoms[j]);
			m_qt_sco.insert(std::map<std::vector<int>, double>::value_type(
				qt_pair, q_lareas[q_suratoms[i]].tc1(t_lareas[t_suratoms[j]]) ) );
		}
	}
	// superpose
	std::pair<int,int> p;
	std::vector<std::pair<int,int> > list_p0;
	std::vector<double> list_score0;
	Position c_pos,g_vec;
	std::vector<Position> list_cpos0,list_gvec0;
	std::map<int, std::vector<int> >::iterator it = smap.begin();
	while( it != smap.end())
	{
		int qid = (*it).first;
		for(int j=0;j < (*it).second.size();j++){
			int tid = (*it).second[j];
			double score = bestsuperpose_larea(q_lareas[qid],t_lareas[tid],c_pos,g_vec);
			if(score >= m_score_threshold){
				p.first = qid;p.second = tid;
				list_p0.push_back(p);
				list_score0.push_back(score);
				list_cpos0.push_back(c_pos);
				list_gvec0.push_back(g_vec);
			}
		}
		++it;
	}

	// sort pair-score list
	std::vector<int> qnum(q_prot.natoms(),0);
	m_list_pair.clear();
	m_list_cpos.clear();
	m_list_gvec.clear();
	m_list_score.clear();
	std::vector<double> list_score1 = list_score0;
	std::vector<bool> chk(list_p0.size(),true);
	std::sort(list_score1.begin(),list_score1.end(),std::greater<double>());
	for(std::vector<double>::iterator it = list_score1.begin(); it != list_score1.end(); ++it){
		for(int i = 0; i < list_p0.size(); i++){
			if(chk[i] && list_score0[i] == (*it)){
			if(qnum[list_p0[i].first] < m_list_pair_maximum){
				m_list_pair.push_back(list_p0[i]);
				m_list_cpos.push_back(list_cpos0[i]);
				m_list_gvec.push_back(list_gvec0[i]);
				m_list_score.push_back(list_score0[i]);
				chk[i] = false;
				qnum[list_p0[i].first] += 1;
			}}
		}
	}
}

//clustering
void AltPS::clustering(const std::vector<LocArea>& q_lareas, const std::vector<LocArea>& t_lareas)
{
	m_list_clusterid.clear();
	m_clustermem_q.clear();
	m_clustermem_t.clear();
	// set first member (single cluster)
	m_clustermem_q.resize(m_list_pair.size());
	m_clustermem_t.resize(m_list_pair.size());
	for(int i = 0; i < m_list_pair.size(); i++) {
		m_clustermem_q[i].push_back(m_list_pair[i].first);
		m_clustermem_t[i].push_back(m_list_pair[i].second);
		m_list_clusterid.push_back(i);
	}
	// clustering main
	for(int i = 0; i < m_list_pair.size()-1; i++){
		for(int j = i+1; j < m_list_pair.size(); j++){
			if(!check_share_member(m_list_clusterid[i],m_list_clusterid[j]))
					if(check_distance(i,j,q_lareas,t_lareas)) cluster2(i,j);
		}
	}
}

//check_share_member
bool AltPS::check_share_member(int clid1, int clid2)
{
	for(int i = 0; i < m_clustermem_q[clid1].size(); i++) {
		for(int j = 0; j < m_clustermem_q[clid2].size(); j++) {
			if(m_clustermem_q[clid1][i] == m_clustermem_q[clid2][j]) return true;
		}
	}
	for(int i = 0; i < m_clustermem_t[clid1].size(); i++) {
		for(int j = 0; j < m_clustermem_t[clid2].size(); j++) {
			if(m_clustermem_t[clid1][i] == m_clustermem_t[clid2][j]) return true;
		}
	}
	return false;
}

//check_cluster_member
bool AltPS::check_cluster_member(std::string qt, int clid, int mid)
{
	for(int i = 0; i < m_list_pair.size(); i++) {
		if(m_list_clusterid[i] == clid){
			int memid;
			if(qt == "q"){memid = m_list_pair[i].first;}
			else{memid = m_list_pair[i].second;}
			if(memid == mid) return true;
		}
	}
	return false;	
}

//check_distance
bool AltPS::check_distance(int a, int b, const std::vector<LocArea>& q_lareas, const std::vector<LocArea>& t_lareas)
{
	LocArea qa = q_lareas[m_list_pair[a].first ];
	LocArea qb = q_lareas[m_list_pair[b].first ];
	LocArea ta = t_lareas[m_list_pair[a].second];
	LocArea tb = t_lareas[m_list_pair[b].second];
	double len_qab1 = m_dmat_q[m_list_pair[a].first][m_list_pair[b].first];
	double len_qab2 = m_list_cpos[a].distance_from(m_list_cpos[b]);
	if( len_qab1 > 10.0 ) return false;
	if( len_qab2 > 10.0 ) return false;
	if(fabs(len_qab1 - len_qab2) > 4.0) return false;
	if(fabs( acos(angle(qa.g_vec, qb.g_vec))
		- acos(angle(m_list_gvec[a],m_list_gvec[b])) ) > 0.523586) return false;// >=PI/6
	LocArea qarea = qa.mix(qb);
	LocArea tarea = ta.mix(tb);
	int qc = qarea.get_atom_number(qb.c_atom);
	int tc = tarea.get_atom_number(tb.c_atom);
	if( bestsuperpose_couple(qarea,tarea,qc,tc) >= 0.70 ) return true;
	return false;
}

//bestsuperpose_couple
double AltPS::bestsuperpose_couple(LocArea& lq, LocArea& lt, int qc, int tc)
{
	Position lqc_pos, lqg_vec;
	// preparation
	std::map<int,int> qid_p;//atom_id -> order of s_atoms
	std::map<int,int> tid_p;
	for(int i = 0;i < lq.nsatoms();i++) qid_p.insert( std::pair<int,int>(lq.s_atoms[i],i) );
	for(int i = 0;i < lt.nsatoms();i++) tid_p.insert( std::pair<int,int>(lt.s_atoms[i],i) );

	std::vector<std::vector<double> > qts(lq.nsatoms());
	for(int i = 0; i< lq.nsatoms(); i++) qts[i].resize(lt.nsatoms(),0.0);
	std::map<std::vector<int>, double>::const_iterator it;
	std::vector<int> qtid(2);
	for(int i = 0;i < lq.nsatoms();i++){
		qtid[0] = lq.s_atoms[i];
		for(int j = 0;j < lt.nsatoms();j++){
			qtid[1] = lt.s_atoms[j];
			it = m_qt_sco.find(qtid);
			if(it != m_qt_sco.end()) qts[ qid_p[qtid[0]] ][ tid_p[qtid[1]] ] = it->second;
		}
	}
	//
	std::vector<std::vector<int> > q_topl;
	for(int i = 0; i< lq.nsatoms(); i++){
		std::vector<int> tlist;
		tlist = set_highscore_target(i,qts);
		q_topl.push_back(tlist);
	}
	//
	// calculate superpose
	std::vector<Position> lq_pos = lq.satom_positions;
	std::vector<Position> lt_pos = lt.satom_positions;
	double dis_th = 1.0;
	double qdis0,qdis1,qdis2;
	double tdis0,tdis1,tdis2;
	double score = -1.0;
	double score_max = -1.0;
	Superposer spp;
	int i = 0;int j = qc;
	for(int k = 1; k< lq.nsatoms(); k++){
		if(k != j){
			qdis1 = m_dmat_q[ lq.s_atoms[i] ][ lq.s_atoms[k] ];
			qdis2 = m_dmat_q[ lq.s_atoms[j] ][ lq.s_atoms[k] ];
			int ti = 0;int tj = tc;
			for(int kk = 0; kk< q_topl[k].size(); kk++){
				int tk = q_topl[k][kk];
				tdis1 = m_dmat_t[ lt.s_atoms[ti] ][ lt.s_atoms[tk] ];
				tdis2 = m_dmat_t[ lt.s_atoms[tj] ][ lt.s_atoms[tk] ];
				if(tj != tk  && fabs(qdis1 - tdis1) <= dis_th && fabs(qdis2 - tdis2) <= dis_th){
					score = superpose_cluster(lq_pos,lt_pos,lq.g_pos,lt.g_vec,i,j,k,ti,tj,tk,qts,spp); 
					if(score > score_max){score_max = score;lq.spp = spp;}
				}
			}
		}
	}
	return score_max;
}

//bestsuperpose_larea
double AltPS::bestsuperpose_larea(LocArea& lq, LocArea& lt, Position& lqc_pos, Position& lqg_vec)
{
	// preparation
	std::map<int,int> qid_p;//atom_id -> order of s_atoms
	std::map<int,int> tid_p;
	for(int i = 0;i < lq.nsatoms();i++) qid_p.insert( std::pair<int,int>(lq.s_atoms[i],i) );
	for(int i = 0;i < lt.nsatoms();i++) tid_p.insert( std::pair<int,int>(lt.s_atoms[i],i) );

	std::vector<std::vector<double> > qts(lq.nsatoms());
	for(int i = 0; i< lq.nsatoms(); i++) qts[i].resize(lt.nsatoms(),0.0);
	std::map<std::vector<int>, double>::const_iterator it;
	std::vector<int> qtid(2);
	for(int i = 0;i < lq.nsatoms();i++){
		qtid[0] = lq.s_atoms[i];
		for(int j = 0;j < lt.nsatoms();j++){
			qtid[1] = lt.s_atoms[j];
			it = m_qt_sco.find(qtid);
			if(it != m_qt_sco.end()) qts[ qid_p[qtid[0]] ][ tid_p[qtid[1]] ] = it->second;
		}
	}
	//
	std::vector<std::vector<int> > q_topl;
	for(int i = 0; i< lq.nsatoms(); i++){
		std::vector<int> tlist;
		tlist = set_highscore_target(i,qts);
		q_topl.push_back(tlist);
	}
	// calculate superpose
	std::vector<Position> lq_pos = lq.satom_positions;
	std::vector<Position> lt_pos = lt.satom_positions;
	double dis_th = 1.0;
	double qdis0,qdis1,qdis2;
	double tdis0,tdis1,tdis2;
	double score = -1.0;
	double score_max = -1.0;
	Superposer spp;
	Position c_pos,g_vec;
	for(int i = 0; i< lq.nsatoms()-2; i++){
	for(int j = i+1; j< lq.nsatoms()-1; j++){
		qdis0 = m_dmat_q[ lq.s_atoms[i] ][ lq.s_atoms[j] ];
	for(int k = j+1; k< lq.nsatoms(); k++){
		qdis1 = m_dmat_q[ lq.s_atoms[i] ][ lq.s_atoms[k] ];
		qdis2 = m_dmat_q[ lq.s_atoms[j] ][ lq.s_atoms[k] ];
	for(int ii = 0; ii< q_topl[i].size(); ii++){int ti = q_topl[i][ii];
	for(int jj = 0; jj< q_topl[j].size(); jj++){int tj = q_topl[j][jj];
		tdis0 = m_dmat_t[ lt.s_atoms[ti] ][ lt.s_atoms[tj] ];
		if(fabs(qdis0 - tdis0) <= dis_th && ti != tj){
	for(int kk = 0; kk< q_topl[k].size(); kk++){int tk = q_topl[k][kk];
		tdis1 = m_dmat_t[ lt.s_atoms[ti] ][ lt.s_atoms[tk] ];
		tdis2 = m_dmat_t[ lt.s_atoms[tj] ][ lt.s_atoms[tk] ];
		if(tj != tk && tk != ti && fabs(qdis1 - tdis1) <= dis_th && fabs(qdis2 - tdis2) <= dis_th){
			score = superpose_larea(lq_pos,lt_pos,lq.g_pos,lt.g_pos,
					i,j,k,ti,tj,tk,qts,c_pos,g_vec,lt.g_vec,spp); 
			if(score > score_max){score_max = score;lqc_pos = c_pos;lqg_vec = g_vec;lq.spp = spp;}
		}
	}}}}
	}}}
	return score_max;
}

//set_highscore_target
std::vector<int> AltPS::set_highscore_target(int i, const std::vector<std::vector<double> >& qts)
{
	std::vector<int> tlist;
	for(int j = 0;j < qts[i].size();j++){
		if (qts[i][j] >= 0.85) tlist.push_back(j);
	}
	return tlist;
}

//set_highscore_target
std::vector<int> AltPS::set_highscore_target(int i, const std::vector<std::vector<double> >& qts, int topnum)
{
	std::vector<int> tlist0;
	std::vector<double> slist0;
	for(int j = 0;j < qts[i].size();j++){
		if (qts[i][j] >= 0.85){
			tlist0.push_back(j);
			slist0.push_back(qts[i][j]);
		}
	}
	int num = std::min(topnum, (int)slist0.size());
	std::vector<int> tlist;
	std::vector<double> slist = slist0;
	std::partial_sort( slist.begin(), slist.begin() + num, slist.end(), std::greater<double>() );
	for(int j = 0; j < tlist0.size(); j++){
		if(slist0[j] >= slist[num-1]) tlist.push_back(tlist0[j]);
		if(tlist.size() == topnum) return tlist;
	}
	return tlist;
}

//superpose_cluster
double AltPS::superpose_cluster(std::vector<Position> a, const std::vector<Position>& b,
	Position ag, const Position& bv, int i, int j, int k,int ti, int tj, int tk,
	const std::vector<std::vector<double> >& qts, Superposer& spp)
{
	double score = 0.0;
	std::vector<Position> pa;
	pa.push_back(a[i]);pa.push_back(a[j]);pa.push_back(a[k]);
	std::vector<Position> pb;
	pb.push_back(b[ti]);pb.push_back(b[tj]);pb.push_back(b[tk]);
	//Superposer spp;
	spp.superpose(pa,pb);
	for (std::vector<Position>::iterator it = a.begin(); it != a.end(); ++it)
		spp.transform(*it);
	spp.transform(ag);
	//calculate score;
	Position av = a[0] - ag;
	if(angle(av,bv) >= 0.0){
		score = simple_alignment(a,b,qts);	
	}else{
		score = 0.0;
	}
	return score;
}

//superpose_larea
double AltPS::superpose_larea(std::vector<Position> a, const std::vector<Position>& b,
	Position ag, const Position& bg, int i, int j, int k,int ti, int tj, int tk,
	const std::vector<std::vector<double> >& qts, Position& cpos, Position& av, const Position& bv, Superposer& spp)
{
	double score = 0.0;
	std::vector<Position> pa;
	pa.push_back(a[i]);pa.push_back(a[j]);pa.push_back(a[k]);
	std::vector<Position> pb;
	pb.push_back(b[ti]);pb.push_back(b[tj]);pb.push_back(b[tk]);
	//Superposer spp;
	spp.superpose(pa,pb);
	for (std::vector<Position>::iterator it = a.begin(); it != a.end(); ++it)
		spp.transform(*it);
	spp.transform(ag);
	cpos = a[0];
	//calculate score;
	av = a[0] - ag;
	if(angle(av,bv) >= 0.0 && cpos.is_within_range(b[0],2.0)){
		score = simple_alignment(a,b,qts);	
	}else{
		score = 0.0;
	}
	return score;
}

//simple_alignment
double AltPS::simple_alignment(std::vector<Position>& a, const std::vector<Position>& b,
    const std::vector<std::vector<double> >& qts)
{
	double score = 0.0;
	for(int i = 0;i < a.size(); i++){
		double smax = 0.0;
		for(int j = 0;j < b.size(); j++){
			if(a[i].is_within_range(b[j], 2.0)){
				if(qts[i][j] > smax) smax = qts[i][j];
			}
		}
		score += smax;
	}
	if(a.size() == 0) return 0.0;
	return score/a.size();
}

//cluster2
void AltPS::cluster2(int a, int b)
{
	int min_id = std::min(m_list_clusterid[a],m_list_clusterid[b]); 	
	int max_id = std::max(m_list_clusterid[a],m_list_clusterid[b]); 	
	for(int i = 0; i < m_list_pair.size(); i++) {
		if(m_list_clusterid[i] == max_id) m_list_clusterid[i] = min_id;
	}
	for(int i = 0; i < m_clustermem_q[max_id].size(); i++)
		m_clustermem_q[min_id].push_back(m_clustermem_q[max_id][i]);
	for(int i = 0; i < m_clustermem_t[max_id].size(); i++)
		m_clustermem_t[min_id].push_back(m_clustermem_t[max_id][i]);
	 m_clustermem_q[max_id].clear();
	 m_clustermem_t[max_id].clear();
}

//analyze_result
void AltPS::analyze_result(const Protein& q_prot, const std::vector<LocArea>& q_lareas,
	const Protein& t_prot, const std::vector<LocArea>& t_lareas)
{
	make_larea_of_each_cluster(q_lareas,t_lareas);
	alignment_cluster_areas(q_prot,t_prot);
}

//make_larea_of_each_cluster
void  AltPS::make_larea_of_each_cluster(const std::vector<LocArea>& q_lareas,
	const std::vector<LocArea>& t_lareas)
{
	m_clusterarea_q.clear();
	m_clusterarea_t.clear();
	for(int clid = 0; clid < m_clustermem_q.size(); clid++){
		std::map<int,Position> map_cl_q;
		for(int i = 0; i < m_clustermem_q[clid].size(); i++){
			int qid = m_clustermem_q[clid][i];
			map_lareas(map_cl_q,q_lareas[qid]);
		}
		std::map<int,Position> map_cl_t;
		for(int i = 0; i < m_clustermem_t[clid].size(); i++){
			int tid = m_clustermem_t[clid][i];
			map_lareas(map_cl_t,t_lareas[tid]);
		}
		m_clusterarea_q.push_back(map_cl_q);
		m_clusterarea_t.push_back(map_cl_t);
	}
}

//map_lareas
void AltPS::map_lareas(std::map<int,Position>& lmap, const LocArea& l)
{
	for(int i = 0; i < l.nsatoms(); i++){
		lmap.insert( std::pair<int, Position>( l.s_atoms[i],l.satom_positions[i]) );	
	}
}

//alignment_cluster_areas
void AltPS::alignment_cluster_areas(const Protein& q_prot,const Protein& t_prot)
{
	std::string basename = m_dirname + m_pairname +"_"+ m_mol1_id +"_"+ m_mol2_id + "_region";
	std::string f_list = basename + ".list";
	std::ofstream fo0(f_list.c_str());
	fo0 << "query  pdbfile [q] = " << m_infile1 << std::endl;
	fo0 << "target pdbfile [t] = " << m_infile2 << std::endl;
	fo0 << "Region_id\tAtom_number_of_detected_area[q]\tAtom_number_of_detected_area[t]\t";
	fo0 << "Alignment_atom_number\tSimilarity-score(S)\tZ-score" << std::endl;
	//
	int outid = 0;
	for(int clid = 0; clid < m_list_clusterid.size(); clid++){
	std::map<int,int> cl_catoms_map;
	for(int i = 0; i < m_list_clusterid.size(); i++){
		if(m_list_clusterid[i] == clid){
			cl_catoms_map.insert(m_list_pair[i]);	
		}
	}
	std::map<int,int> cl_catoms_map_nr = non_redundant(cl_catoms_map);
	if(cl_catoms_map_nr.size() >= m_cluster_minimum){
		Superposer* spp;
		spp = calc_cluster_superposer(clid, cl_catoms_map_nr);

		LocArea cl_area_q = get_cluster_larea(clid,"q");
		LocArea cl_area_t = get_cluster_larea(clid,"t");
		double score_sum = superpose_cluster_area(clid,cl_area_q,cl_area_t,(*spp),q_prot,t_prot);
		int min_num = std::min(cl_area_q.nsatoms(), cl_area_t.nsatoms());
		double score;if(min_num > 0){score = score_sum / min_num;}else{score = 0.0;}
		double zscore = calc_zscore(min_num,score);

		//---  print_out  [./dataset/results/] ---
		if(zscore >= m_zscore && min_num >= m_cluster_natom){
			outid += 1;
			fo0 << outid << "\t" << cl_area_q.nsatoms() << "\t" << cl_area_t.nsatoms() << "\t" 
				<< m_alignment_atom_list.size() << "\t" << std::setprecision(2) << score << "\t" 
				<< std::setprecision(3) << zscore <<  std::endl;
			//[m_pairname_m_mol1_id_m_mol2_id_region_id_alignment_atoms.list]
			std::string f_alignment_list = basename + "_" + i_to_s(outid) +"_alignment_atoms.list";
			std::ofstream fo1(f_alignment_list.c_str());
			std::map<int, int>::iterator it = m_alignment_atom_list.begin();
			while( it != m_alignment_atom_list.end() )
			{
				int qid = (*it).first;
				int tid = (*it).second;
				int q_serial = q_prot.atom(qid)->number;
				int t_serial = t_prot.atom(tid)->number;
				const Atom* qatom = q_prot.atom(qid);
				const Atom* tatom = t_prot.atom(tid);
				int qclass = q_prot.atom(qid)->pc_class;
				int tclass = t_prot.atom(tid)->pc_class;
				const Protein::Residue qres = q_prot.residue(m_resnumbers_q[qid]);
				const Protein::Residue tres = t_prot.residue(m_resnumbers_t[tid]);
				fo1 << q_serial << "\t" << qatom->name<< "\t" << qres.name 
					<< "\t" << qres.number << "\t" << patty_atname(qclass) << "\t:\t"
					<< t_serial << "\t" << tatom->name<< "\t" << tres.name
					<< "\t" << tres.number << "\t" << patty_atname(tclass) << std::endl;
				++it;
			}
			//[m_pairname_m_mol1_id_m_mol2_id_region_id_area_q.pdb]
			std::string f_pdb1 = basename + "_" + i_to_s(outid) +"_area_q.pdb";
			cl_area_q.spp = (*spp);
			cl_area_q.output_superpose_pdb(q_prot,f_pdb1,false);
			//[m_pairname_m_mol1_id_m_mol2_id_region_id_area_t.pdb]
			std::string f_pdb2 = basename + "_" + i_to_s(outid) +"_area_t.pdb";
			cl_area_t.output_pdb(t_prot,f_pdb2,false);
			//[m_pairname_m_mol1_id_m_mol2_id_region_id_protein_q.pdb]
			std::string f_pdb3 = basename + "_" + i_to_s(outid) +"_protein_q.pdb";
			output_superpose_pdb(q_prot,f_pdb3,*spp);
			//[m_pairname_m_mol1_id_m_mol2_id_region_id_protein_t.pdb]
			std::string f_pdb4 = basename + "_" + i_to_s(outid) +"_protein_t.pdb";
			output_pdb(t_prot,f_pdb4);
			//[m_pairname_m_mol1_id_m_mol2_id_region_id_spp.dat]
			std::string f_pdb5 = basename + "_" + i_to_s(outid) +"_spp.dat";
			std::ofstream fo5(f_pdb5.c_str());
			(*spp).print(fo5);
		}
	}}
	std::cout << outid << " regions are detected. [Z-score >= " << m_zscore << ", Size >= " << m_cluster_natom << "]" << std::endl;
	std::cout << "region list file: " << f_list << std::endl;
	std::cout << std::endl;
}

//non_redundant
std::map<int,int> AltPS::non_redundant(const std::map<int,int>& cl_catoms_map)
{
	std::map<int, int> cl_catoms_map_nr;
	std::map<int, int>::const_iterator it = cl_catoms_map.begin();
	while( it != cl_catoms_map.end() )
	{
		cl_catoms_map_nr.insert(std::pair<int,int>((*it).second, (*it).first));
		++it;
	}
	return cl_catoms_map_nr;
}

//calc_cluster_superposer
Superposer* AltPS::calc_cluster_superposer(int clid, const std::map<int,int>& cl_catoms_map)
{
	std::vector<Position> pa;
	std::vector<Position> pb;
	std::map<int, int>::const_iterator it = cl_catoms_map.begin();
	while( it != cl_catoms_map.end() )
	{
		pa.push_back(m_clusterarea_q[clid][(*it).second]);
		pb.push_back(m_clusterarea_t[clid][(*it).first]);
		++it;
	}
	static Superposer spp;
	spp.superpose(pa,pb);
	return &spp;
}

//get_cluster_larea
LocArea AltPS::get_cluster_larea(int clid, std::string qt)
{
	LocArea cl_larea;
	if(qt == "q"){
		std::map<int,Position> clusterarea =  m_clusterarea_q[clid];
		std::map<int,Position>::iterator it = clusterarea.begin();
		while( it != clusterarea.end() )
		{
			int atom_id = (*it).first;
			cl_larea.s_atoms.push_back(atom_id);
			cl_larea.resn_of_satoms.push_back(m_resnumbers_q[atom_id]);
			cl_larea.satom_positions.push_back((*it).second);
			++it;
		}
	}else{
		std::map<int,Position> clusterarea =  m_clusterarea_t[clid];
		std::map<int,Position>::iterator it = clusterarea.begin();
		while( it != clusterarea.end() )
		{
			int atom_id = (*it).first;
			cl_larea.s_atoms.push_back(atom_id);
			cl_larea.resn_of_satoms.push_back(m_resnumbers_t[atom_id]);
			cl_larea.satom_positions.push_back((*it).second);
			++it;
		}
	}
	return cl_larea;
}

//superpose_cluster_area
double AltPS::superpose_cluster_area(int clid, LocArea& cl_area_q, LocArea& cl_area_t, const Superposer& spp,
	const Protein& q_prot,const Protein& t_prot)
{
	m_alignment_atom_list.clear();

	for(int i = 0; i < cl_area_q.satom_positions.size(); i++){
		spp.transform(cl_area_q.satom_positions[i]);
	}

	std::vector<std::vector<double> > len(cl_area_q.nsatoms());//distance between atom_i of area_q and atom_j of area_t
	for(int i = 0; i < cl_area_q.nsatoms(); ++i){
		for(int j = 0; j < cl_area_t.nsatoms(); ++j){
			len[i].push_back(cl_area_q.satom_positions[i].distance_from(cl_area_t.satom_positions[j]));
		}
	}
	//
	std::vector<int> q_pair(cl_area_q.nsatoms(),-1);
	std::vector<int> t_pair(cl_area_t.nsatoms(),-1);
	int near = 1;
	while(near == 1){
		near  = 0;
		std::vector<int> q__t(cl_area_q.nsatoms(),-1);
		for(int i = 0; i < cl_area_q.nsatoms(); ++i){
		if(q_pair[i] == -1){
			double dlentmp = 2.5;
			for(int j = 0; j < cl_area_t.nsatoms(); ++j){
			if(t_pair[j] == -1){
				if(len[i][j] <= dlentmp){
					if( patty_match( (q_prot.atom(cl_area_q.s_atoms[i]))->pc_class,
									 (t_prot.atom(cl_area_t.s_atoms[j]))->pc_class ))  
					{q__t[i] = j;dlentmp = len[i][j];near = 1;}
				}
			}}
		}}
		std::vector<int> t__q(cl_area_t.nsatoms(),-1);
		for(int i = 0; i < cl_area_t.nsatoms(); ++i){
		if(t_pair[i] == -1){
			double dlentmp = 2.5;
			for(int j = 0; j < cl_area_q.nsatoms(); ++j){
			if(q_pair[j] == -1){
				if(len[j][i] <= dlentmp){
					if( patty_match( (q_prot.atom(cl_area_q.s_atoms[j]))->pc_class,
									 (t_prot.atom(cl_area_t.s_atoms[i]))->pc_class  ) )  
					{t__q[i] = j;dlentmp = len[j][i];near = 1;}
				}
			}}
		}}
		for(int i = 0; i < cl_area_q.nsatoms(); ++i){
			if(q_pair[i] == -1){
				int jj = q__t[i];
				if(jj >= 0 && t_pair[jj] == -1 && t__q[jj] == i){
					q_pair[i] = jj;
					t_pair[jj] = i;
				}
			}
		}
	}

	//score
	double score = 0.0;
	for(int i = 0; i < cl_area_q.nsatoms(); ++i){
		if(q_pair[i] >= 0){
			m_alignment_atom_list.insert(std::pair<int,int>(cl_area_q.s_atoms[i], cl_area_t.s_atoms[q_pair[i]]));
			std::vector<int> qt_pair;
			qt_pair.push_back(cl_area_q.s_atoms[i]);
			qt_pair.push_back(cl_area_t.s_atoms[q_pair[i]]);
			score += m_qt_sco[qt_pair];
		}
	}
	return score;
}

//calc_zscore
double AltPS::calc_zscore(int n, double score)
{
	double ave = 0.733 / std::pow(n,0.30588);
	double sd = 0.686 / std::pow(n,0.596);
	return (score - ave) / sd;
}
