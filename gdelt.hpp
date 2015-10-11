#include <string>

template<typename T>
struct DiffClass {
  T value_;
  double distance(const DiffClass& other) {
    return fabs(value_ - other.value_);
  }
};

template<typename T>
struct CategoryClass {
  T value_;
  double distance(const CategoryClass& other) {
    return value_ == other.value_;
  }
};

template<typename T>
std::istream& operator>>(std::istream& in, DiffClass<T>& v) {
  return in >> v.value_;
}

template<typename T>
std::istream& operator>>(std::istream& in, CategoryClass<T>& v) {
  return in >> v.value_;
}

typedef DiffClass<int> DaysSinceEpoch;
typedef DiffClass<int> EventCode;
typedef DiffClass<int> QuadClass;
typedef DiffClass<double> GoldsteinScale;
typedef CategoryClass<std::string> Actor1Name;
typedef CategoryClass<std::string> Actor2Name;
typedef DiffClass<bool> IsRootEvent;
typedef DiffClass<int> NumMentions;
typedef DiffClass<int> NumSources;
typedef DiffClass<int> NumArticles;
typedef DiffClass<double> AvgTone;
typedef CategoryClass<std::string> Actor1Geo;
typedef CategoryClass<std::string> Actor2Geo;
typedef CategoryClass<std::string> BaseURL;

typedef std::tuple<DaysSinceEpoch, EventCode, QuadClass,
                   GoldsteinScale, Actor1Name, Actor2Name,
                   IsRootEvent, NumMentions, NumSources, NumArticles,
                   AvgTone, Actor1Geo, Actor2Geo, BaseURL> GDELTMini;
