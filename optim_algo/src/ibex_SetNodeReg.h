//============================================================================
//                                  I B E X                                   
// File        : ibex_SetNode.h
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : 13 juil. 2014
//============================================================================

#ifndef __IBEX_SET_NODE_REG_H__
#define __IBEX_SET_NODE_REG_H__

#include "ibex.h"
#include "ibex_NodeType.h"
#include <unistd.h>

using namespace std;
namespace ibex {

static long compteur=0;
/**
 * \brief Set node.
 */
class SetNodeReg {
public: 
	SetNodeReg * right;
	SetNodeReg * left;
	bool cflag;
	NodeType status;
	unsigned int var;
	
	/**
	 * \brief Create two sons to a node.
	 */
	 void cut(const IntervalVector& box);
	 void cut(const IntervalVector& box,unsigned int cutvar);
	 void changeStat(const NodeType& statp, const NodeType& statn);
	 SetNodeReg * getRight();
	 SetNodeReg * getLeft();
	 IntervalVector left_box(const IntervalVector& nodebox) const;

	IntervalVector right_box(const IntervalVector& nodebox) const;
	void gatherFrom(int type);
	bool is_leaf() const;

private:
	friend class SetIntervalReg;
	//friend std::ostream& operator<<(std::ostream& os, const SetIntervalReg& set);
	/**
	 * \brief Callback for "visit_leaves"
	 */
	typedef void (*leaf_func) (const IntervalVector&, BoolInterval);

	/**
	 * \brief Creates a node of given status
	 */
	SetNodeReg(NodeType status);

	/**
	 * \brief Creates a node of given status and father
	 */
	 SetNodeReg(NodeType status,SetNodeReg *father);

	 /**
	 * \brief Creates a node of given status and cut variable, use only to copy a node
	 */
	 SetNodeReg(NodeType status, unsigned int var);

	 /**
	 * \brief Creates a node of given status father and cut variable, use only to cut a node
	 */
	 SetNodeReg(NodeType status,SetNodeReg *father, unsigned int var);

	 SetNodeReg(NodeType status, SetNodeReg *father,  unsigned int var, double pt);

	/**
	 * \brief Delete this.
	 */
	~SetNodeReg();

	/**
	 * \brief True iff this node is a leaf.
	 */
	

	/**
	 * \brief Return a copy of the node.
	 */
	SetNodeReg * copy( SetNodeReg * nodeFather) const;

	

	/**
	 * \brief Intersection or Union between a box of NodeType val and a regular setInterval.
	 * Box may not respect the regularity of the regular setInterval.
	 * The part of the set that intersect nodebox only is modify.
	 */
	void operator_ir(const IntervalVector& box,const IntervalVector& nodebox, NodeType val, bool op, double eps);

	/**
	 * \brief Overloaded function, apply valout to the set outside of the box and valin to the inside, thus modify the whole set.
	 */
	void operator_ir(const IntervalVector& box,const IntervalVector& subbox, NodeType valin,NodeType valout, bool op,double eps);

	/**
	 * \brief Intersection between two SetNodeReg.
	 */
	void oper(SetNodeReg* other,bool op);

	/**
	 * \brief Intersection between SetNodeReg and a status.
	 */
	void inter(NodeType x_status);

	/**
	 * \brief Union between SetNodeReg and a status.
	 */
	void _union(NodeType x);

	/**
	 * \brief Change value of branch if right and left got the same value IN, OUT or UNK.
	 */
	void gather();

	void gather(int type);

	void gatherFrom();

	


	/**
	 * \brief Change value of leaves of status statp to a new status statn
	 */
	 

	/**
	 * \brief Contrat i-set w.r.t separator sep, split leaves until boxin and boxout returned by sep are disjoints, 
	 * and then call operator_ir to contract on them
	 */	
	void cleave(const IntervalVector& box, Sep& sep, const double eps);

	void cleave(const IntervalVector& box, Ctc& Ctcin, Ctc& Ctcout, const double eps);

	void findNeighbor(const IntervalVector& box,const IntervalVector& nbox,vector<SetNodeReg*> * neigh) ;

	void findNeighbor(const IntervalVector& box,const IntervalVector& nbox,vector<SetNodeReg*> * neigh, vector<IntervalVector> * neighbox) ;

	void findBox(const IntervalVector& bbox,IntervalVector * box);

	void compare(const IntervalVector& root, const IntervalVector& box,const int& type, bool * flag, unsigned int op);

	void getLeaf(vector<SetNodeReg*> * vecLeaf);
	
	void getLeaf(vector<SetNodeReg*> * vecLeaf, vector<IntervalVector> * vecBox, const IntervalVector& curbox);

	/**
	 * \brief Visit the leaves of a tree with a callback "func"
	 */
	void visit_leaves(leaf_func func, const IntervalVector& nodebox) const;

	/**
	 * \brief Display the structure on output stream "os"
	 */
	void print(std::ostream& os, const IntervalVector& nodebox, int shift) const;

	




	/**
	 * \brief The status of the node
	 */
	SetNodeReg * father;
	
	double pt;
	
	
};

int Interset(IntervalVector x, IntervalVector y);
//std::ostream& operator<<(std::ostream& os, const SetIntervalIr& set);

bool is_point(IntervalVector v); 
} // namespace ibex

#endif // __IBEX_SET_NODE_REG_H__
