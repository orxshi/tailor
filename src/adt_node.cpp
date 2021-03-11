#include "adt.h"

namespace Tailor
{
    void Node::rotate(double ang, int axis, const Vector3& rot_point)
    {
        if (p_ != nullptr)
            p_->rotate(ang, axis, rot_point);

        if (left_ != nullptr)
            left_->rotate(ang, axis, rot_point);

        if (right_ != nullptr)
            right_->rotate(ang, axis, rot_point);
    }

    void Node::move(const Vector3& v)
    {
        if (p_ != nullptr)
            p_->move(v);

        if (left_ != nullptr)
            left_->move(v);

        if (right_ != nullptr)
            right_->move(v);
    }

    double Node::c(int i) const
    {
        assert(i >= 0);
        assert(i < TAILOR_ADT_VAR);
        return c_[i];
    }

    double Node::d(int i) const
    {
        assert(i >= 0);
        assert(i < TAILOR_ADT_VAR);
        return d_[i];
    }

    void Node::construct_id_address_map(std::map<int, Node*>& map)
    {
        if (p_ != nullptr)
        {
            map.insert(std::make_pair(p_->idx(), this));
        }

        if (left_ != nullptr)
            left_->construct_id_address_map(map);
        if (right_ != nullptr)
            right_->construct_id_address_map(map);
    }

    Node::Node():
        p_(nullptr),
        key_(0.),
        left_(nullptr),
        right_(nullptr),
        level_(0),
        dim_(0)
    {
    }

    Node::Node(const Node& other): c_(other.c_), level_(other.level_), d_(other.d_), key_(other.key_), left_(nullptr), dim_(other.dim_), p_(nullptr), right_(nullptr)
    {
        //std::cout << "copy cons of node" << std::endl;
        if (other.p_ != nullptr)
        {
            //delete p_;
            p_ = new ADTPoint(*other.p_);
        }

        if (other.left_ != nullptr)
        {
            //delete left_;
            left_ = new Node(*other.left_);
        }

        if (other.right_ != nullptr)
        {
            //delete right_;
            right_ = new Node(*other.right_);
        }
        //std::cout << "exit copy cons of node" << std::endl;
    }

    Node::Node(unsigned int level, const std::vector<double>& c, const std::vector<double>& d): left_(nullptr), p_(nullptr), right_(nullptr), c_(c), level_(level), d_(d)
    {
        dim_ = level_ % TAILOR_ADT_VAR;
        key_ = 0.5 * (c[dim_] + d[dim_]);
    }

    Node::Node(unsigned int level, const std::vector<double>& c, const std::vector<double>& d, const ADTPoint& p): Node(level, c, d)
    {
        p_ = new ADTPoint(p);
    }

    Node::~Node()
    {    
        destroy_children();
    }

    void Node::destroy_children()
    {
        if (p_ != nullptr)
            delete p_;

        if (left_ != nullptr)
            left_->destroy_children();

        if (right_ != nullptr)
            right_->destroy_children();
    }

    unsigned int Node::level() const
    {
        return level_;
    }

    double Node::key() const
    {
        return key_;
    }

    Node* Node::insert(const ADTPoint& point)
    {
        Node* node = nullptr;

        if (point.dim(dim_) < key_)
        {
            bool added = insert_left(point);
            if (added)
                return left_;
            node = left_->insert(point);
        }
        else
        {
            bool added = insert_right(point);
            if (added)
                return right_;
            node = right_->insert(point);
        }

        return node;
    }

    void Node::add_point(const ADTPoint& point)
    {
        p_ = new ADTPoint(point);
        //isEmpty_ = false;
    }

    bool Node::insert_left(const ADTPoint& point)
    {
        if (left_ == nullptr)
        {
            std::vector<double> newd = d_;

            for (int d=0; d<TAILOR_ADT_VAR; ++d)
            {
                if (d == dim_)
                {
                    newd[d] = key_;
                }
            }

            left_ = new Node(level_ + 1, c_, newd);
        }

        if (left_->p_ == nullptr)
        {
            left_->add_point(point);
            return true;
        }

        return false;
    }

    bool Node::insert_right(const ADTPoint& point)
    {
        if (right_ == nullptr)
        {
            std::vector<double> newc = c_;

            for (int d=0; d<TAILOR_ADT_VAR; ++d)
            {
                if (d == dim_)
                {
                    newc[d] = key_;
                }
            }

            right_ = new Node(level_ + 1, newc, d_);
        }

        if (right_->p_ == nullptr)
        {
            right_->add_point(point);
            return true;
        }

        return false;
    }

    bool Node::doRegionOverlapCheckAll(const ADTPoint& targetPoint) const
    {
        assert(false);
        /*//std::cout << "checck all" << std::endl;
        assert(c_.size() == TAILOR_ADT_VAR);
        assert(d_.size() == TAILOR_ADT_VAR);
        assert(targetPoint.dim().size() == TAILOR_ADT_VAR);

        for (int d=0; d<TAILOR_ADT_DIM; ++d)
        {
            double dis1 = c_[2*d] - targetPoint.dim(2*d+1);
            double dis2 = d_[2*d+1] - targetPoint.dim(2*d);
            //if (c_[2*d] > targetPoint.dim(2*d+1))
            if (dis1 > TAILOR_ZERO)
                return false;
            //if (d_[2*d+1] < targetPoint.dim(2*d))
            if (dis2 < TAILOR_ZERO)
                return false;
        }
        //std::cout << "end checck all" << std::endl;

        return true;*/

        //Point p(targetPoint.dim(0), targetPoint.dim(2));
        //return aabb_.do_intersect(p);
    }

    bool Node::doRegionOverlap(const ADTPoint& targetPoint, bool verbose) const
    {
        if ((dim_ % 2) == 0)
        {
            double dis1 = c_[dim_] - targetPoint.dim(dim_+1);
            if (verbose)
            {
                std::cout << "doRegionOverlap dim_ | c_[dim_] | targetPoint.dim(dim_+1) | dis1 | overlap " << dim_ << " " << c_[dim_] << " " << targetPoint.dim(dim_+1) << " " << dis1 << " " << !(dis1 > TAILOR_ZERO) << std::endl;
            }
            //if (c_[dim_] > targetPoint.dim(dim_+1))
            if (dis1 > TAILOR_ZERO)
            //if (dis1 > 0.)
                return false;
        }
        else
        {
            //double dis1 = d_[dim_] - targetPoint.dim(dim_-1);
            double dis1 = targetPoint.dim(dim_-1) - d_[dim_];
            if (verbose)
            {
                std::cout << "doRegionOverlap dim_ | d_[dim_] | targetPoint.dim(dim_-1) | dis1 | overlap " << dim_ << " " << d_[dim_] << " " << targetPoint.dim(dim_-1) << " " << dis1 << " " << !(dis1 < TAILOR_ZERO) << std::endl;
            }
            //if (d_[dim_] < targetPoint.dim(dim_-1))
            //if (dis1 < TAILOR_ZERO)
            if (dis1 > TAILOR_ZERO)
            //if (dis1 < 0.)
                return false;
        }

        return true;
    }

    bool Node::overlap(const ADTPoint& targetPoint, bool verbose) const
    {
        //bool verbose = false;

        //if (targetPoint.idx() == 537)
        //verbose = true;
        //else
        //verbose = false;

        //if (p_ == nullptr) return false;
        //if (!doCubesOverlap(targetPoint)) return false;
        //if (!compareFunction(targetPoint)) return false;

        if (p_ != nullptr)
        {
            if (verbose)
                std::cout << "p is not null" << std::endl;
            if (p_->overlap(targetPoint, verbose)) return true;
        }
        else
        {
            if (verbose)
                std::cout << "p is null" << std::endl;
        }

        return false;
    }

    void Node::searchChildren(const ADTPoint& targetPoint, std::vector<const Node*>& searchStack, bool verbose) const
    {
        if (left_ != nullptr)
        {
            if (verbose)
            {
                std::cout << "checking left" << std::endl;
            }
            if (left_->doRegionOverlap(targetPoint, verbose))
            {
                searchStack.push_back(left_);
            }
        }

        if (right_ != nullptr)
        {
            if (verbose)
            {
                std::cout << "checking right" << std::endl;
            }
            if (right_->doRegionOverlap(targetPoint, verbose))
            {
                searchStack.push_back(right_);
            }
        }
    }

    const ADTPoint* Node::p() const
    {
        return p_;
    }

    void Node::remove_point()
    {
        assert(p_ != nullptr);
        delete p_;
        p_ = nullptr;
    }
}
