#ifndef QUADTREE_H
#define QUADTREE_H

#include "core.h"
#include "aabb.h"

#define QT_DIM 2

struct QTPoint
{
    int id;
    Point p;
    QTPoint(int id, Point p): id(id), p(p) {};
};

struct QTNode
{
    AABB aabb;
    std::array<AABB, 4> children_aabb;
    QTPoint p;
    std::array<QTNode*, 4> quadrant;

    QTNode(QTPoint p, AABB aabb);
};

struct Quadtree
{
    QTNode* root;
    std::vector<QTNode*> searchStack;
    std::vector<int> ids;
    std::vector<QTPoint> points;

    Quadtree(std::vector<QTPoint>& points);

    void destroy_tree (QTNode *&leaf);
    bool insert (QTPoint& point);
    bool insert_(QTPoint& point, QTNode* const& parent, QTNode*& child, int iquad);
    void searchChildren (QTNode* const& node, const QTPoint& targetPoint);
    void search (QTNode* const& node, const QTPoint& targetPoint);
    void search (const QTPoint& targetPoint);
    void build();    
};

#endif

