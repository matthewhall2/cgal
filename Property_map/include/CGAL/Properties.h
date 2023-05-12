#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <CGAL/assertions.h>

#include <map>

// todo: maybe this could be avoided
#include <boost/property_map/property_map.hpp>

namespace CGAL::Properties {

template <typename Index>
class Property_array_base {
public:

  Property_array_base() = default;

  Property_array_base(const Property_array_base<Index>& rhs) = delete;

  virtual ~Property_array_base() = default;

  // Declare virtual functions here, for things which need to be done within the Property container
  // todo: maybe these should be private, and made available using friend

  virtual std::shared_ptr<Property_array_base<Index>> clone(const std::vector<bool>& active_indices) = 0;

  virtual void copy(const Property_array_base<Index>& other) = 0;

  virtual void append(const Property_array_base<Index>& other) = 0;

  virtual void reserve(std::size_t n) = 0;

  virtual void swap(Index a, Index b) = 0;

  virtual void reset(Index i) = 0;

};

/*!
 * \brief Indexed storage for arbitrary types
 *
 * @tparam T
 */
template <typename Index, typename T>
class Property_array : public Property_array_base<Index> {

  std::vector<T> m_data;
  const std::vector<bool>& m_active_indices;
  T m_default_value;

public:

  using value_type = T;
  using reference = typename std::vector<T>::reference;
  using const_reference = typename std::vector<T>::const_reference;

  Property_array(const std::vector<bool>& active_indices, const T& default_value) :
    m_data(), m_active_indices(active_indices), m_default_value(default_value) {

    m_data.reserve(active_indices.capacity());
    m_data.resize(active_indices.size(), m_default_value);
  }

  virtual std::shared_ptr<Property_array_base<Index>> clone(const std::vector<bool>& active_indices) override {
    auto new_array = std::make_shared<Property_array<Index, T>>(active_indices, m_default_value);
    new_array->m_data = m_data;
    return new_array;
  }

  virtual void copy(const Property_array_base<Index>& other_base) {
    auto& other = dynamic_cast<const Property_array<Index, T>&>(other_base);
    m_data = other.m_data;
    CGAL_precondition(m_active_indices.size() == m_data.size());
  }

  virtual void append(const Property_array_base<Index>& other_base) override {
    auto& other = dynamic_cast<const Property_array<Index, T>&>(other_base);
    CGAL_precondition(m_data.size() + other.m_data.size() == m_active_indices.size());
    m_data.insert(m_data.end(), other.m_data.begin(), other.m_data.end());
  }

  virtual void reserve(std::size_t n) override {
    CGAL_precondition(m_active_indices.size() == n);
    m_data.resize(n, m_default_value);
  };

  virtual void swap(Index a, Index b) override {
    // todo: maybe cast to index, instead of casting index to size?
    CGAL_precondition(std::size_t(a) < m_data.size() && std::size_t(b) < m_data.size());
    std::iter_swap(m_data.begin() + a, m_data.begin() + b);
  };

  virtual void reset(Index i) override {
    CGAL_precondition(std::size_t(i) < m_data.size());
    m_data[std::size_t(i)] = m_default_value;
  };

  std::size_t capacity() const { return m_data.size(); }

public:

  const_reference operator[](Index i) const {
    CGAL_precondition(std::size_t(i) < m_data.size());
    return m_data[std::size_t(i)];
  }

  reference operator[](Index i) {
    CGAL_precondition(std::size_t(i) < m_data.size());
    return m_data[std::size_t(i)];
  }

public:

  bool operator==(const Property_array<Index, T>& other) const {
    return &other == this;
  }

  bool operator!=(const Property_array<Index, T>& other) const { return !operator==(other); }

};


template <typename Index, typename T>
class Property_array_handle {

  Property_array<Index, T>& m_array;

public:

  // Necessary for use as a boost::property_type
  using key_type = Index;
  using value_type = T;
  using reference = typename std::vector<T>::reference;
  using const_reference = typename std::vector<T>::const_reference;
  using category = boost::lvalue_property_map_tag;

  Property_array_handle(Property_array<Index, T>& array) : m_array(array) {}

  const_reference operator[](Index i) const { return m_array[i]; }

  reference operator[](Index i) { return m_array[i]; }

  bool operator==(const Property_array<Index, T>& other) const { return &other.m_array == m_array; }

  bool operator!=(const Property_array<Index, T>& other) const { return !operator==(other); }

  inline friend reference get(Property_array_handle<Index, T> p, const Index& i) { return p[i]; }

  inline friend void put(Property_array_handle<Index, T> p, const Index& i, const T& v) { p[i] = v; }

};

template <typename Index = std::size_t>
class Property_container {

  std::map<std::string, std::shared_ptr<Property_array_base<Index>>> m_property_arrays{};
  std::vector<bool> m_active_indices{};

public:

  Property_container() = default;

  Property_container(const Property_container<Index>& other) {
    m_active_indices = other.m_active_indices;
    for (auto [name, array]: other.m_property_arrays) {
      // todo: this could probably be made faster using emplace_hint
      m_property_arrays.emplace(
        name,
        array->clone(m_active_indices)
      );
    }
  }

  Property_container(Property_container<Index>&& other) { *this = std::move(other); }

  // todo: maybe this could be implemented in terms of the move assignment operator?
  Property_container<Index>& operator=(const Property_container<Index>& other) {
    m_active_indices = other.m_active_indices;
    for (auto [name, other_array]: other.m_property_arrays) {

      // If this container has a property by the same name
      auto it = m_property_arrays.find(name);
      if (it != m_property_arrays.end()) {
        auto [_, this_array] = *it;

        // No naming collisions with different types allowed
        CGAL_precondition(typeid(*this_array) == typeid(*other_array));

        // Copy the data from the other array
        this_array->copy(*other_array);

      } else {
        // Adds the new property
        m_property_arrays.emplace(name, other_array->clone(m_active_indices));
      }
    }
    return *this;
  }

  Property_container<Index>& operator=(Property_container<Index>&& other) {
    m_active_indices = std::move(other.m_active_indices);
    for (auto [name, other_array]: other.m_property_arrays) {

      // If this container has a property by the same name
      auto it = m_property_arrays.find(name);
      if (it != m_property_arrays.end()) {
        auto [_, this_array] = *it;

        // No naming collisions with different types allowed
        CGAL_precondition(typeid(*this_array) == typeid(*other_array));

        // Copy the data from the other array
        this_array->copy(*other_array);

      } else {
        // Adds the new property
        m_property_arrays.emplace(name, other_array->clone(m_active_indices));
      }
    }

    // The moved-from property map should retain all of its properties, but contain 0 elements
    other.reserve(0);
    return *this;
  }

  template <typename T>
  std::pair<std::reference_wrapper<Property_array<Index, T>>, bool>
  get_or_add_property(const std::string& name, const T default_value = T()) {
    auto [it, created] = m_property_arrays.emplace(
      name,
      std::make_shared<Property_array<Index, T>>(
        m_active_indices,
        default_value
      )
    );
    auto [key, array] = *it;
    auto& typed_array = dynamic_cast<Property_array<Index, T>&>(*array);
    return {{typed_array}, created};
  }

  template <typename T>
  Property_array<Index, T>& add_property(const std::string& name, const T default_value = T()) {
    // todo: I'm not settled on the naming, but it's really convenient to have a function like this
    auto [array, created] = get_or_add_property(name, default_value);
    CGAL_precondition(created);
    return array.get();
  }

  template <typename T>
  const Property_array<Index, T>& get_property(const std::string& name) const {
    CGAL_precondition(m_property_arrays.count(name) != 0);
    return dynamic_cast<const Property_array<Index, T>&>(*m_property_arrays.at(name));
  }

  template <typename T>
  Property_array<Index, T>& get_property(const std::string& name) {
    CGAL_precondition(m_property_arrays.count(name) != 0);
    return dynamic_cast<Property_array<Index, T>&>(*m_property_arrays.at(name));
  }

  /*!
   * Removes a property array from the container
   *
   * @param name
   * @return True if a container with this name existed, false otherwise
   */
  bool remove_property(const std::string& name) { return m_property_arrays.erase(name) == 1; }

  std::size_t num_properties() { return m_property_arrays.size(); }

public:

  void reserve(std::size_t n) {
    m_active_indices.resize(n);
    for (auto [name, array]: m_property_arrays)
      array->reserve(n);
  }

  std::size_t size() const { return std::count(m_active_indices.begin(), m_active_indices.end(), true); }

  std::size_t capacity() const { return m_active_indices.size(); }

  Index emplace_back() {

    // Expand the storage and return the last element
    reserve(capacity() + 1);
    m_active_indices.back() = true;
    Index first_new_index{capacity() - 1};
    reset(first_new_index);
    return first_new_index;
  }

  Index emplace() {

    // If there are empty slots, return the index of one of them and mark it as full
    auto first_unused = std::find_if(m_active_indices.begin(), m_active_indices.end(), [](bool used) { return !used; });
    if (first_unused != m_active_indices.end()) {
      *first_unused = true;
      auto index = Index(std::distance(m_active_indices.begin(), first_unused));
      reset(index);
      return index;
    }

    return emplace_back();
  }

  Index emplace_group_back(std::size_t n) {

    // Expand the storage and return the start of the new region
    reserve(capacity() + n);
    for (auto it = m_active_indices.end() - n; it < m_active_indices.end(); ++it)
      *it = true;
    return Index(capacity() - n);
  }

  Index emplace_group(std::size_t n) {

    auto search_start = m_active_indices.begin();
    while (search_start != m_active_indices.end()) {

      // Find the first unused cell
      auto unused_begin = std::find_if(
        search_start, m_active_indices.end(),
        [](bool used) { return !used; }
      );

      // Determine if the group fits
      auto unused_end = std::find_if(
        unused_begin, std::min(unused_begin + n, m_active_indices.end()),
        [](bool used) { return used; }
      );

      // If the discovered range was large enough
      if (std::distance(unused_begin, unused_end) >= n) {

        // Mark the indices as used, and reset the properties of each of them
        // todo: it would be better to provide a function to set a range
        for (auto it = unused_begin; it < unused_end; ++it) {
          *it = true;
          reset(Index(std::distance(m_active_indices.begin(), it)));
        }

        // Return the first index of the range
        return Index(std::distance(m_active_indices.begin(), unused_begin));
      }

      // If we didn't find a large enough region, continue our search after the end
      search_start = unused_end;
    }

    // If no empty regions were found, expand the storage
    return emplace_group_back(n);
  }

  void swap(Index a, Index b) {
    for (auto [name, array]: m_property_arrays)
      array->swap(a, b);
  }

  void reset(Index i) {
    for (auto [name, array]: m_property_arrays)
      array->reset(i);
  }

  void erase(Index i) {
    m_active_indices[i] = false;
    for (auto [name, array]: m_property_arrays)
      array->reset(i);
  }

  bool is_erased(Index i) const {
    return !m_active_indices[i];
  }

  // todo: I'd prefer to eliminate this, if possible
  void mark_active(Index i) {
    return m_active_indices[i] = true;
  }

  /*!
   * Adds the elements of the other container to this container for each property which is present in this container.
   *
   * Gaps are preserved, and all elements of the other container are guaranteed
   * to appear after the elements of this container.
   * Properties in this container which don't appear in the other container are extended with default values.
   * Properties in the other container which don't appear in this one are not included.
   * todo: merge() would be useful as well, but could break contiguous regions in the other container
   *
   * @param other
   */
  void append(const Property_container<Index>& other) {
    // todo

    m_active_indices.insert(m_active_indices.end(), other.m_active_indices.begin(), other.m_active_indices.end());
    for (auto [name, array]: m_property_arrays) {
      auto it = other.m_property_arrays.find(name);
      if (it != other.m_property_arrays.end())
        array->append(*it->second);
      else
        array->reserve(m_active_indices.size());
    }
  }
};

}

#endif //ORTHTREE_TESTS_PROPERTIES_H
