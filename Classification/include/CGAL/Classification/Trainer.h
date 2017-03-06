// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_TRAINER_H
#define CGAL_CLASSIFICATION_TRAINER_H

#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Label.h>
#include <CGAL/Classifier.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#include <tbb/mutex.h>
#endif // CGAL_LINKED_WITH_TBB

//#define CGAL_CLASSTRAINING_VERBOSE
#if defined(CGAL_CLASSTRAINING_VERBOSE)
#define CGAL_CLASSTRAINING_CERR std::cerr
#else
#define CGAL_CLASSTRAINING_CERR std::ostream(0)
#endif

namespace CGAL {

namespace Classification {

/*!
\ingroup PkgClassification

\brief Training algorithm to set up weights and effects of the
features used for classification.

Each label must have ben given a small set of user-defined inliers to
provide the training algorithm with a ground truth (see `set_inlier()`
and `set_inliers()`).

This methods estimates the set of feature weights and of
[effects](@ref Classification::Feature::Effect) that make the
classifier succeed in correctly classifying the sets of inliers given
by the user. These parameters are directly modified within the
`Classification::Feature_base` and `Classification::Label` objects.

\tparam ItemRange model of `ConstRange`. Its iterator type is
`RandomAccessIterator`.

\tparam ItemMap model of `ReadablePropertyMap` whose key
type is the value type of the iterator of `ItemRange` and value type is
the type of the items that are classified.

\tparam ConcurrencyTag enables sequential versus parallel
algorithm. Possible values are `Parallel_tag` (default value is CGAL
is linked with TBB) or `Sequential_tag` (default value otherwise).
*/
template <typename ItemRange, typename ItemMap,
#if defined(DOXYGEN_RUNNING)
          typename ConcurrencyTag>
#elif defined(CGAL_LINKED_WITH_TBB)
          typename ConcurrencyTag = CGAL::Parallel_tag>
#else
          typename ConcurrencyTag = CGAL::Sequential_tag>
#endif
class Trainer
{

public:
  typedef CGAL::Classifier<ItemRange, ItemMap, ConcurrencyTag>  Classifier;
  typedef typename Classification::Label_handle                 Label_handle;
  typedef typename Classification::Feature_handle               Feature_handle;

private:

#ifdef CGAL_LINKED_WITH_TBB
  class Compute_worst_score_and_confidence
  {
    std::vector<std::size_t>& m_training_set;
    Classifier* m_classifier;
    std::size_t m_label;
    double& m_confidence;
    std::size_t& m_nb_okay;
    tbb::mutex& m_mutex;
    
  public:

    Compute_worst_score_and_confidence (std::vector<std::size_t>& training_set,
                                        Classifier* classifier,
                                        std::size_t label,
                                        double& confidence,
                                        std::size_t& nb_okay,
                                        tbb::mutex& mutex)
      : m_training_set (training_set)
      , m_classifier (classifier)
      , m_label (label)
      , m_confidence (confidence)
      , m_nb_okay (nb_okay)
      , m_mutex (mutex)
    { }
          
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t k = r.begin(); k != r.end(); ++ k)
        {
          std::vector<std::pair<double, std::size_t> > values;
      
          for(std::size_t l = 0; l < m_classifier->number_of_labels(); ++ l)
            values.push_back (std::make_pair (m_classifier->energy_of (m_classifier->label(l),
                                                                       m_training_set[k]),
                                              l));
          std::sort (values.begin(), values.end());

          if (values[0].second == m_label)
            {
              m_mutex.lock();
              m_confidence += values[1].first - values[0].first;
              ++ m_nb_okay;
              m_mutex.unlock();
            }
        }
    }

  };
#endif // CGAL_LINKED_WITH_TBB

  
  Classifier* m_classifier;
  std::vector<std::vector<std::size_t> > m_training_sets;
  std::vector<double> m_precision;
  std::vector<double> m_recall;
  std::vector<double> m_iou; // intersection over union
  double m_accuracy;
  double m_mean_iou;
  double m_mean_f1;

public:

  /// \name Constructor
  /// @{
  
  Trainer (Classifier& classifier)
    : m_classifier (&classifier)
    , m_training_sets (classifier.number_of_labels())
  {
  }

  /// @}


  /// \name Inliers
  /// @{
  
  /*!
    \brief Adds the item at position `index` as an inlier of
    `label` for the training algorithm.

    \note This inlier is only used for training. There is no guarantee
    that the item at position `index` will be classified as `label`
    after calling `run()`, `run_with_local_smoothing()` or
    `run_with_graphcut()`.

    \return `true` if the inlier was correctly added, `false`
    otherwise (if `label` was not found).
  */
  bool set_inlier (Label_handle label, std::size_t index)
  {
    std::size_t label_idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_classifier->number_of_labels(); ++ i)
      if (m_classifier->label(i) == label)
        {
          label_idx = i;
          break;
        }
    if (label_idx == (std::size_t)(-1))
      return false;

    if (label_idx >= m_training_sets.size())
      m_training_sets.resize (label_idx + 1);

    m_training_sets[label_idx].push_back (index);

    return true;
  }

  /*!

    \brief Adds the items at positions `indices` as inliers of
    `label` for the training algorithm.

    \note These inliers are only used for training. There is no
    guarantee that the items at positions `indices` will be classified
    as `label` after calling `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()`.

    \tparam IndexRange range of `std::size_t`, model of `ConstRange`.
  */
  template <class IndexRange>
  bool set_inliers (Label_handle label,
                    IndexRange indices)
  {
    std::size_t label_idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_classifier->number_of_labels(); ++ i)
      if (m_classifier->label(i) == label)
        {
          label_idx = i;
          break;
        }
    if (label_idx == (std::size_t)(-1))
      return false;

    if (label_idx >= m_training_sets.size())
      m_training_sets.resize (label_idx + 1);

    std::copy (indices.begin(), indices.end(), std::back_inserter (m_training_sets[label_idx]));

    return true;
  }
  
  /*!
    \brief Resets inlier sets used for training.
  */
  void reset_inlier_sets()
  {
    std::vector<std::vector<std::size_t > >().swap (m_training_sets);
  }

  
  /// @}


  /// \name Training
  /// @{

  /*!
    \brief Runs the training algorithm.

    All the `Classification::Label` and `Classification::Feature`
    necessary for classification should have been added before running
    this function. After training, the user can call `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()` to compute
    the classification using the estimated parameters.

    \param nb_tests number of tests to perform. Higher values may
    provide the user with better results at the cost of a higher
    computation time. Using a value of at least 10 times the number of
    features is advised.

    \return minimum ratio (over all labels) of provided
    ground truth items correctly classified using the best
    configuration found.
  */
  
  double train (std::size_t nb_tests = 300)
  {
    if (m_training_sets.empty())
      return 0.;

    for (std::size_t i = 0; i < m_classifier->number_of_labels(); ++ i)
      if (m_training_sets.size() <= i || m_training_sets[i].empty())
        std::cerr << "WARNING: \"" << m_classifier->label(i)->name() << "\" doesn't have a training set." << std::endl;

    std::vector<double> best_weights (m_classifier->number_of_features(), 1.);

    struct Feature_training
    {
      bool skipped;
      double wmin;
      double wmax;
      double factor;
    };
    std::vector<Feature_training> feature_train;
    std::size_t nb_trials = 100;
    double wmin = 1e-5, wmax = 1e5;
    double factor = std::pow (wmax/wmin, 1. / (double)nb_trials);
    std::size_t feature_used = 0;

    std::vector<std::pair<std::size_t, std::size_t> > sorted_features;
    
    for (std::size_t j = 0; j < m_classifier->number_of_features(); ++ j)
      {
        Feature_handle feature = m_classifier->feature(j);
        best_weights[j] = feature->weight();

        std::size_t nb_useful = 0;
        double min = (std::numeric_limits<double>::max)();
        double max = -(std::numeric_limits<double>::max)();

        feature->set_weight(wmin);
        for (std::size_t i = 0; i < 100; ++ i)
          {
            estimate_feature_effect(feature);
            if (feature_useful(feature))
              {
                CGAL_CLASSTRAINING_CERR << "#";
                nb_useful ++;
                min = (std::min) (min, feature->weight());
                max = (std::max) (max, feature->weight());
              }
            else
              CGAL_CLASSTRAINING_CERR << "-";
            feature->set_weight(factor * feature->weight());
          }
        CGAL_CLASSTRAINING_CERR << std::endl;
        CGAL_CLASSTRAINING_CERR << feature->name() << " useful in "
                                << nb_useful << "% of the cases, in interval [ "
                                << min << " ; " << max << " ]" << std::endl;
        feature_train.push_back (Feature_training());
        feature_train.back().skipped = false;
        feature_train.back().wmin = min / factor;
        feature_train.back().wmax = max * factor;
        if (nb_useful < 2)
          {
            feature_train.back().skipped = true;
            feature->set_weight(0.);
            best_weights[j] = feature->weight();
          }
        else if (best_weights[j] == 1.)
          {
            feature->set_weight(0.5 * (feature_train.back().wmin + feature_train.back().wmax));
            best_weights[j] = feature->weight();
            sorted_features.push_back (std::make_pair (-nb_useful, j));
            ++ feature_used;
          }
        else
          {
            feature->set_weight(best_weights[j]);
            sorted_features.push_back (std::make_pair (-nb_useful, j));
            ++ feature_used;
          }
        estimate_feature_effect(feature);
      }

    std::size_t nb_trials_per_feature = 1 + (std::size_t)(nb_tests / (double)(feature_used));
    CGAL_CLASSIFICATION_CERR << "Trials = " << nb_tests << ", features = " << feature_used
              << ", trials per feature = " << nb_trials_per_feature << std::endl;
    for (std::size_t i = 0; i < feature_train.size(); ++ i)
      if (!(feature_train[i].skipped))
        feature_train[i].factor = std::pow (feature_train[i].wmax / feature_train[i].wmin,
                                        1. / (double)nb_trials_per_feature);
    
    
    double best_score = 0.;
    double best_confidence = 0.;
    boost::tie (best_confidence, best_score)
      = compute_worst_confidence_and_score (0., 0.);
    
    CGAL_CLASSIFICATION_CERR << "TRAINING GLOBALLY: Best score evolution: " << std::endl;

    CGAL_CLASSIFICATION_CERR << 100. * best_score << "% (found at initialization)" << std::endl;

    std::size_t current_feature_changed = 0;
    CGAL::Timer teff, tconf, tscore;

    std::sort (sorted_features.begin(), sorted_features.end());
    for (std::size_t i = 0; i < sorted_features.size(); ++ i)
      {
        current_feature_changed = sorted_features[i].second;
        // while (feature_train[current_feature_changed].skipped)
        //   {
        //     ++ current_feature_changed;
        //     if (current_feature_changed == m_classifier->number_of_features())
        //       current_feature_changed = 0;
        //   }

        std::size_t nb_used = 0;
        for (std::size_t j = 0; j < m_classifier->number_of_features(); ++ j)
          {
            if (j == current_feature_changed)
              continue;
            
            m_classifier->feature(j)->set_weight(best_weights[j]);
            teff.start();
            estimate_feature_effect(m_classifier->feature(j));
            teff.stop();
            if (feature_useful(m_classifier->feature(j)))
              nb_used ++;
          }
        Feature_handle current_feature = m_classifier->feature(current_feature_changed);
        const Feature_training& tr = feature_train[current_feature_changed];
        
        current_feature->set_weight(tr.wmin);
        for (std::size_t j = 0; j < nb_trials_per_feature; ++ j)
          {
            teff.start();
            estimate_feature_effect(current_feature);
            teff.stop();

            tconf.start();
            double worst_confidence = 0., worst_score = 0.;
            boost::tie (worst_confidence, worst_score)
              = compute_worst_confidence_and_score (best_confidence, best_score);
            tconf.stop();

            if (worst_score > best_score
                && worst_confidence > best_confidence)
              {
                best_score = worst_score;
                best_confidence = worst_confidence;
                CGAL_CLASSIFICATION_CERR << 100. * best_score << "% (found at iteration "
                          << (i * nb_trials_per_feature) + j << "/" << nb_tests << ", "
                          << nb_used + (feature_useful(current_feature) ? 1 : 0)
                          << "/" << m_classifier->number_of_features() << " feature(s) used)" << std::endl;
                for (std::size_t k = 0; k < m_classifier->number_of_features(); ++ k)
                  {
                    Feature_handle feature = m_classifier->feature(k);
                    best_weights[k] = feature->weight();
                  }
              }
            current_feature->set_weight(current_feature->weight() * tr.factor);
          }

        ++ current_feature_changed;
      }

    std::cerr << "Estimation of effects = " << teff.time() << std::endl
              << "Confidence computation = " << tconf.time() << std::endl
              << "Score computation = " << tscore.time() << std::endl;
    
    for (std::size_t i = 0; i < best_weights.size(); ++ i)
      {
        Feature_handle feature = m_classifier->feature(i);
        feature->set_weight(best_weights[i]);
      }

    estimate_features_effects();
    
    CGAL_CLASSIFICATION_CERR << std::endl << "Best score found is at least " << 100. * best_score
              << "% of correct classification" << std::endl;

    std::size_t nb_removed = 0;
    for (std::size_t i = 0; i < best_weights.size(); ++ i)
      {
        Feature_handle feature = m_classifier->feature(i);
        CGAL_CLASSTRAINING_CERR << "FEATURE " << feature->name() << ": " << best_weights[i] << std::endl;
        feature->set_weight(best_weights[i]);

        Classification::Feature::Effect side = m_classifier->label(0)->feature_effect(feature);
        bool to_remove = true;
        for (std::size_t j = 0; j < m_classifier->number_of_labels(); ++ j)
          {
            Label_handle clabel = m_classifier->label(j);
            if (clabel->feature_effect(feature) == Classification::Feature::FAVORING)
              CGAL_CLASSTRAINING_CERR << " * Favored for ";
            else if (clabel->feature_effect(feature) == Classification::Feature::PENALIZING)
              CGAL_CLASSTRAINING_CERR << " * Penalized for ";
            else
              CGAL_CLASSTRAINING_CERR << " * Neutral for ";
            if (clabel->feature_effect(feature) != side)
              to_remove = false;
            CGAL_CLASSTRAINING_CERR << clabel->name() << std::endl;
          }
        if (to_remove)
          {
            CGAL_CLASSTRAINING_CERR << "   -> Useless! Should be removed" << std::endl;
            ++ nb_removed;
          }
      }
    CGAL_CLASSIFICATION_CERR << nb_removed
              << " feature(s) out of " << m_classifier->number_of_features() << " are useless" << std::endl;

    compute_precision_recall();
    return best_score;
  }

  /// @}


  /// \name Evaluation
  /// @{

  /*!

    \brief Returns the precision of the training for the given label.

    Precision is the number of true positives divided by the sum of
    the true positives and the false positives.

  */
  double precision (Label_handle label) const
  {
    std::size_t label_idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_classifier->number_of_labels(); ++ i)
      if (m_classifier->label(i) == label)
        {
          label_idx = i;
          break;
        }
    if (label_idx == (std::size_t)(-1))
      return 0.;

    return m_precision[label_idx];
  }

  /*!

    \brief Returns the recall of the training for the given label.

    Recall is the number of true positives divided by the sum of
    the true positives and the false negatives.

  */
  double recall (Label_handle label) const
  {
    std::size_t label_idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_classifier->number_of_labels(); ++ i)
      if (m_classifier->label(i) == label)
        {
          label_idx = i;
          break;
        }
    if (label_idx == (std::size_t)(-1))
      return 0.;

    return m_recall[label_idx];
  }

  /*!

    \brief Returns the \f$F_1\f$ score of the training for the given label.

    \f$F_1\f$ score is the harmonic mean of `precision()` and `recall()`:

    \f[
    F_1 = 2 \times \frac{precision \times recall}{precision + recall}
    \f]

  */
  double f1_score (Label_handle label) const
  {
    std::size_t label_idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_classifier->number_of_labels(); ++ i)
      if (m_classifier->label(i) == label)
        {
          label_idx = i;
          break;
        }
    if (label_idx == (std::size_t)(-1))
      return 0.;
    
    return 2. * (m_precision[label_idx] * m_recall[label_idx])
      / (m_precision[label_idx] + m_recall[label_idx]);
  }

  /*!
    \brief Returns the intersection over union of the training for the
    given label.

    Intersection over union is the number of true positives divided by
    the sum of the true positives, of the false positives and of the
    false negatives.
  */
  double intersection_over_union (Label_handle label) const
  {
    std::size_t label_idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_classifier->number_of_labels(); ++ i)
      if (m_classifier->label(i) == label)
        {
          label_idx = i;
          break;
        }
    if (label_idx == (std::size_t)(-1))
      return 0.;
    
    return m_iou[label_idx];
  }

  /*!
    \brief Returns the accuracy of the training.

    Accuracy is the total number of true positives divided by the
    total number of provided inliers.
  */
  double accuracy() const { return m_accuracy; }
  
  /*!
    \brief Returns the mean \f$F_1\f$ score of the training over all
    labels (see `f1_score()`).
  */
  double mean_f1_score() const { return m_mean_f1; }
  
  /*!
    \brief Returns the mean intersection over union of the training
    over all labels (see `intersection_over_union()`).
  */
  double mean_intersection_over_union() const { return m_mean_iou; }
  
  /// @}

  /// \cond SKIP_IN_MANUAL
  Label_handle training_label_of (std::size_t index) const
  {
    for (std::size_t i = 0; i < m_training_sets.size(); ++ i)
      if (std::find (m_training_sets[i].begin(), m_training_sets[i].end(), index) != m_training_sets[i].end())
        return m_classifier->label(i);
    return Label_handle();
  }

  /// \endcond

private:

  void estimate_features_effects()
  {
    for (std::size_t i = 0; i < m_classifier->number_of_features(); ++ i)
      estimate_feature_effect (m_classifier->feature(i));
  }

  void estimate_feature_effect (Feature_handle feature)
  {
    std::vector<double> mean (m_classifier->number_of_labels(), 0.);
                                  
    for (std::size_t j = 0; j < m_classifier->number_of_labels(); ++ j)
      {
        for (std::size_t k = 0; k < m_training_sets[j].size(); ++ k)
          {
            double val = feature->normalized(m_training_sets[j][k]);
            mean[j] += val;
          }
        mean[j] /= m_training_sets[j].size();
      }

    std::vector<double> sd (m_classifier->number_of_labels(), 0.);
        
    for (std::size_t j = 0; j < m_classifier->number_of_labels(); ++ j)
      {
        Label_handle clabel = m_classifier->label(j);
            
        for (std::size_t k = 0; k < m_training_sets[j].size(); ++ k)
          {
            double val = feature->normalized(m_training_sets[j][k]);
            sd[j] += (val - mean[j]) * (val - mean[j]);
          }
        sd[j] = std::sqrt (sd[j] / m_training_sets[j].size());
        if (mean[j] - sd[j] > 0.75)
          clabel->set_feature_effect (feature, Classification::Feature::FAVORING);
        else if (mean[j] + sd[j] < 0.25)
          clabel->set_feature_effect (feature, Classification::Feature::PENALIZING);
        else
          clabel->set_feature_effect (feature, Classification::Feature::NEUTRAL);
      }
  }

  std::pair<double, double> compute_worst_confidence_and_score (double lower_conf, double lower_score)
  {
    double worst_confidence = (std::numeric_limits<double>::max)();
    double worst_score = (std::numeric_limits<double>::max)();
    
    for (std::size_t j = 0; j < m_classifier->number_of_labels(); ++ j)
      {
        double confidence = 0.;
        std::size_t nb_okay = 0;

#ifndef CGAL_LINKED_WITH_TBB
        CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                                   "Parallel_tag is enabled but TBB is unavailable.");
#else
        if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
          {
            tbb::mutex mutex;
            Compute_worst_score_and_confidence f(m_training_sets[j], m_classifier, j, confidence, nb_okay, mutex);
            tbb::parallel_for(tbb::blocked_range<size_t>(0, m_training_sets[j].size ()), f);
          }
        else
#endif
          {
            for (std::size_t k = 0; k < m_training_sets[j].size(); ++ k)
              {
                std::vector<std::pair<double, std::size_t> > values;
      
                for(std::size_t l = 0; l < m_classifier->number_of_labels(); ++ l)
                  values.push_back (std::make_pair (m_classifier->energy_of (m_classifier->label(l),
                                                                             m_training_sets[j][k]),
                                                    l));
                std::sort (values.begin(), values.end());

                if (values[0].second == j)
                  {
                    confidence += values[1].first - values[0].first;
                    ++ nb_okay;
                  }
              }
          }
        
        double score = nb_okay / (double)(m_training_sets[j].size());
        confidence /= (double)(m_training_sets[j].size() * m_classifier->number_of_features());

        if (confidence < worst_confidence)
          worst_confidence = confidence;
        if (score < worst_score)
          worst_score = score;
        
        if (worst_confidence < lower_conf || worst_score < lower_score)
          return std::make_pair (worst_confidence, worst_score);
      }
    return std::make_pair (worst_confidence, worst_score);
  }

  bool feature_useful (Feature_handle feature)
  {
    Classification::Feature::Effect side = m_classifier->label(0)->feature_effect(feature);
    for (std::size_t k = 1; k < m_classifier->number_of_labels(); ++ k)
      if (m_classifier->label(k)->feature_effect(feature) != side)
        return true;
    return false;
  }

  
  void compute_precision_recall ()
  {
    std::vector<std::size_t> true_positives (m_classifier->number_of_labels());
    std::vector<std::size_t> false_positives (m_classifier->number_of_labels());
    std::vector<std::size_t> false_negatives (m_classifier->number_of_labels());

    std::size_t sum_true_positives = 0;
    std::size_t total = 0;
    
    for (std::size_t j = 0; j < m_classifier->number_of_labels(); ++ j)
      {
        for (std::size_t k = 0; k < m_training_sets[j].size(); ++ k)
          {
            std::size_t nb_class_best=0; 
            double val_class_best = (std::numeric_limits<double>::max)();
      
            for(std::size_t l = 0; l < m_classifier->number_of_labels(); ++ l)
              {
                double value = m_classifier->energy_of (m_classifier->label(l),
                                                        m_training_sets[j][k]);
          
                if(val_class_best > value)
                  {
                    val_class_best = value;
                    nb_class_best = l;
                  }
              }
            ++ total;
            if (nb_class_best == j)
              {
                ++ true_positives[j];
                ++ sum_true_positives;
              }
            else
              {
                ++ false_negatives[j];
                for(std::size_t l = 0; l < m_classifier->number_of_labels(); ++ l)
                  if (nb_class_best == l)
                    ++ false_positives[l];
              }
          }
      }


    m_precision.clear();
    m_recall.clear();

    m_mean_iou = 0.;
    m_mean_f1 = 0.;
    
    for (std::size_t j = 0; j < m_classifier->number_of_labels(); ++ j)
      {
        m_precision.push_back (true_positives[j] / double(true_positives[j] + false_positives[j]));
        m_recall.push_back (true_positives[j] / double(true_positives[j] + false_negatives[j]));
        m_iou.push_back (true_positives[j] / double(true_positives[j] + false_positives[j] + false_negatives[j]));

        m_mean_iou += m_iou.back();
        m_mean_f1 += 2. * (m_precision.back() * m_recall.back())
          / (m_precision.back() + m_recall.back());
      }

    m_mean_iou /= m_classifier->number_of_labels();
    m_mean_f1 /= m_classifier->number_of_labels();
    m_accuracy = sum_true_positives / double(total);

  }

  

};


} // namespace Classification
  
} // namespace CGAL

#endif // CGAL_CLASSIFICATION_TRAINER_H
