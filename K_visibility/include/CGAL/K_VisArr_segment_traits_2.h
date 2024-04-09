#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

template <typename Kernel_ = CGAL::Exact_predicates_exact_constructions_kernel>
class K_VisArr_segment;

template <typename Kernel_ = CGAL::Exact_predicates_exact_constructions_kernel>
class K_VisArr_segment_traits : public CGAL::Arr_segment_traits_2<Kernel_> {
    
    friend class K_VisArr_segment<Kernel_>;
  //  typedef typename Arr_segment_2<Kernel_> K_VisArr_segment;
  //  typedef typename K_VisArr_segment Arr_segment;
public:
    typedef CGAL::Arr_segment_traits_2<Kernel_> Superclass;

    typedef typename Kernel_::Point_2        Point_2;
    typedef  K_VisArr_segment<Kernel_>           X_monotone_curve_2;
    typedef  K_VisArr_segment<Kernel_>           Curve_2;
};

template <typename Kernel_>
class K_VisArr_segment : public K_VisArr_segment_traits<Kernel_>::_Segment_cached_2 {

    typedef typename K_VisArr_segment_traits<Kernel_>::_Segment_cached_2 Base;
    typedef typename Kernel_::Segment_2                               Segment_2;
    typedef typename Kernel_::Point_2                                 Point_2;
    typedef typename Kernel_::Line_2                                  Line_2;
    private:
        int id_;
public:


    K_VisArr_segment(const typename Segment_2& seg) : Base(seg) {
        this.id_ = seg.id();
    }

    int id() const { return id_; }

    int& id() { return id_; }

        
};