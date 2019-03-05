#ifndef CBA_ALTPS_H
#define CBA_ALTPS_H

#include <functional>
#include "bio_sequence.h"
#include "app.h"
#include "locarea.h"

namespace Cba {

class AltPS : public App {
public:
	AltPS();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	double m_score_smap;
	double m_score_threshold;
	int m_list_pair_maximum;
	int m_cluster_minimum;
	int m_cluster_natom;
	double m_zscore;
	std::string m_dirname;
	std::string m_infile1;
	std::string m_infile2;
	std::string m_pairname;
	std::string m_mol1_id;
	std::string m_mol2_id;
	std::string m_mol1_chainid;
	std::string m_mol2_chainid;
	std::vector<std::vector<double> > m_dmat_q;
	std::vector<std::vector<double> > m_dmat_t;
	std::map<std::vector<int>, double> m_qt_sco;
	std::vector<std::pair<int,int> > m_list_pair;
	std::vector<double> m_list_score;
	std::vector<Position> m_list_cpos;
	std::vector<Position> m_list_gvec;
	std::vector<int> m_list_clusterid;
	std::vector<std::vector<int> > m_clustermem_q;
	std::vector<std::vector<int> > m_clustermem_t;
	std::vector<std::map<int,Position> > m_clusterarea_q;
	std::vector<std::map<int,Position> > m_clusterarea_t;
	std::vector<int> m_resnumbers_q;
	std::vector<int> m_resnumbers_t;
	std::map<int,int> m_alignment_atom_list;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(AltPS* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		AltPS* m_obj;
	};
	friend class SetParam;

	int check_names();
	std::vector<std::vector<double> > calc_locarea(Protein& prot, std::vector<LocArea>& lareas, std::vector<int>& suratoms);
	void compare_with_pdb(Protein& prot, std::vector<LocArea>& lareas, std::vector<int>& suratoms);
	void search_area( 
		const Protein& q_prot, std::vector<LocArea>& q_lareas, std::vector<int>& q_suratoms,
		const Protein& t_prot, std::vector<LocArea>& t_lareas, std::vector<int>& t_suratoms);
	void pairname(const std::string f1,const std::string f2);
	void clustering(const std::vector<LocArea>& q_lareas, const std::vector<LocArea>& t_lareas);
	bool check_cluster_member(std::string qt, int clid, int mid);
	bool check_share_member(int clid1, int clid2);
	bool check_distance(int i, int j, const std::vector<LocArea>& q_lareas, const std::vector<LocArea>& t_lareas);
	void cluster2(int i, int j);
	void analyze_result(
		const Protein& q_prot, const std::vector<LocArea>& q_lareas,
		const Protein& t_prot, const std::vector<LocArea>& t_lareas);
	void make_larea_of_each_cluster(const std::vector<LocArea>& q_lareas,
		const std::vector<LocArea>& t_lareas);
	void map_lareas(std::map<int,Position>& lmap, const LocArea& l);
	void alignment_cluster_areas(const Protein& q_prot,const Protein& t_prot);
	std::map<int,int> non_redundant(const std::map<int,int>& cl_catoms_map);
	Superposer* calc_cluster_superposer(int clid, const std::map<int,int>& cl_catoms_map);
	LocArea get_cluster_larea(int clid, std::string qt);
	double superpose_cluster_area(int clid, LocArea& cl_area_q, LocArea& cl_area_t, const Superposer& spp,
		const Protein& q_prot,const Protein& t_prot);
	double bestsuperpose_couple(LocArea& lq, LocArea& lt, int qc, int tc);
	double superpose_cluster(std::vector<Position> lq_pos, const std::vector<Position>& tq_pos,
	    Position ag, const Position& bv, int i, int j, int k, int ti, int tj, int tk,
	    const std::vector<std::vector<double> >& qts, Superposer& spp);
	double bestsuperpose_larea(LocArea& lq, LocArea& lt, Position& c_pos, Position& g_vec);
	double superpose_larea(std::vector<Position> lq_pos, const std::vector<Position>& tq_pos,
	    Position ag, const Position& bg, int i, int j, int k, int ti, int tj, int tk,
	    const std::vector<std::vector<double> >& qts, Position& cpos, Position& av, const Position& bv, Superposer& spp);
	std::vector<int> set_highscore_target(int i, const std::vector<std::vector<double> >& qts);
	std::vector<int> set_highscore_target(int i, const std::vector<std::vector<double> >& qts, int topnum);
	double simple_alignment(std::vector<Position>& a, const std::vector<Position>& b,
	    const std::vector<std::vector<double> >& qts);
	double calc_zscore(int num, double score);
	std::vector<std::vector<double> > calc_dmat(Protein& prot);
};





}
#endif
