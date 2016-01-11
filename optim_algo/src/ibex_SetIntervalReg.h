//============================================================================
//                                  I B E X                                   
// File        : ibex_Set.h
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : 13 juil. 2014
//============================================================================

#ifndef __IBEX_SET_REG_H__
#define __IBEX_SET_REG_H__

//#include "ibex.h"
#include "ibex_SetIntervalReg.h"
#include "ibex_SetNodeReg.h"

namespace ibex {

/**
 * \defgroup iset Set Interval
 */

/**
 * \ingroup iset
 * \brief Set Interval
 *#include "ibex_SetIntervalReg.h"
#include "ibex_CellHeap.h"
#include "ibex.h"
 * See "Solving set-valued constraint satisfaction problem", Luc Jaulin, Computing. Volume 94, Issue 2, Page 297-311.
 *
 */
class SetIntervalReg{
public:
	double eps;
	SetNodeReg* root; 
	IntervalVector bounding_box; // not sure it is really necessary

      ~SetIntervalReg() {
		delete root;
	}

	/**
	 * \brief Creates a set interval from a simple box
	 *
	 * If x is the bounding box in argument, the set interval is either
	 * [empty, x] if inner==false or [x,x] if inner==false.
	 */
	SetIntervalReg(const IntervalVector& bounding_box, double eps, bool inner=true) ;

	/**
	 * \brief Creates a copy of an existing SetIntervalReg
	 *
	 */
	SetIntervalReg(const SetIntervalReg& set);

	SetIntervalReg(const char *filename);


	/**
	 * \brief Loads a set from a data file.
	 *
	 * \see #save().
	 */
	//SetIntervalReg(const char* filename);



	/**
	 * \brief i-Set Intersection
	 *
	 * In Jaulin's terminology, this operator is the "i-set extension of the intersection".
	 *
	 * If [x] designates this i-set and [y] the i-set in argument, then this will be replace by
	 *  { x \cap y, x\in[x] and y\in[y] }.
	 */
      void operator&=(const SetIntervalReg& set);
      void operator&=(const IntervalVector& set);

	/**
	 * \brief i-Set Intersection
	 *
	 * Intersection of a SetIntervalReg with a box of interior value valin and exterior value valout. 
	 * Computation is more efficient than create a SetIntervalReg contracted on the box, and then
	 * compute intersection of the two SetIntervalReg. 
	 */
      void interBox(const IntervalVector& box,NodeType valin=__IBEX_IN__,NodeType valout=__IBEX_OUT__);


	/**
	 * \brief i-Set Union
	 *
	 * In Jaulin's terminology, this operator is the "i-set extension of the union".
	 *
	 * If [x] designates this i-set and [y] the i-set in argument, then this will be replace by
	 *  { x \cup y, x\in[x] and y\in[y] }.
	 */
      void operator|=(const SetIntervalReg& set);
      void operator|=(const IntervalVector& set);

	/**
	 * \brief i-Set Intersection
	 *
	 * Union of a SetIntervalReg with a box of interior value valin and exterior value valout. 
	 * Computation is more efficient than create a SetIntervalReg contracted on the box, and then
	 * compute intersection of the two SetIntervalReg. 
	 */
      void unionBox(const IntervalVector& box,NodeType valin=__IBEX_IN__,NodeType valout=__IBEX_OUT__);


	/**
	 * \brief True if this i-set is empty
	 *
	 * \warning: an empty i-set is different from a i-set containing (and possibly only containing) the empty set.
	 */
      bool is_empty() const;

      bool is_reg() const;

	/**
	 * \brief Contrat i-set w.r.t a separator sep
	 */
      void contract(Sep& sep);

      void contract(IntervalVector * box,int type);

      void contract(Ctc& Ctcin, Ctc& Ctcout);

      void findNeighbor(const IntervalVector& box, vector<SetNodeReg*> * neigh);

      void findNeighbor(const IntervalVector& box, vector<SetNodeReg*> * neigh,vector<IntervalVector> * neighbox);

      void findNeighbor(SetNodeReg* node, vector<SetNodeReg*> * neigh );

      void findNeighbor(SetNodeReg* node, vector<SetNodeReg*> * neigh, vector<IntervalVector> * neighbox );

      IntervalVector findNodeBox(SetNodeReg * node);

      bool contains(const IntervalVector& box,const int& type) const;

	/**
	 * \brief Check if box intersects the subset of SetInterval of status type;
	 */
      bool intersects(const IntervalVector& box, const int& type) const;

	/**
	 * \brief Check if box is disjoint the subset of SetInterval of status type;
	 */
      bool is_disjoint(const IntervalVector& box,const int& type) const;

      vector<SetNodeReg*> getLeaf() const;

      void getLeaf(vector<SetNodeReg*> * node_vect,vector<IntervalVector> * node_box);

      void gather();

      void gather(int type);

	/**
	 * \brief Serialize the set and save it into a file
	 */
	void save(const char* filename);

      void visit_leaves(SetNodeReg::leaf_func func) const;

	/**
	 * \brief Distance of the point "pt" wrt the set (if inside is true)
	 * of the complementary of the set (if inside is false).
	 */
    //  double dist(const Vector& pt, bool inside) const;

	

	

protected:

	/**
	 * \brief Load the set from a file
	 */
	void load(const char* filename);

	friend std::ostream& operator<<(std::ostream& os, const SetIntervalReg& set);
	
	// NULL means no existing set (warning: different from empty set!)

	


};

std::ostream& operator<<(std::ostream& os, const SetIntervalReg& set);




inline void SetIntervalReg::operator&=(const IntervalVector& set) { interBox(set);};

inline void SetIntervalReg::operator|=(const IntervalVector& set) { unionBox(set);};



} // namespace ibex

#endif // __IBEX_SET_REG_H__
