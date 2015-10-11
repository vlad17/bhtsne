#include <string>

template<typename T>
struct DiffClass {
  T value_;
  double distance(const DiffClass& other) {
    return fabs(value_ - other.value_);
  }
};

struct CategoryClass {
  int64_t value_; // hash of a string or other category
  double distance(const CategoryClass& other) {
    if (value_ == 0 || other.value_ == 0) return 1;
    return value_ != other.value_;
  }
};

template<typename T>
std::istream& operator>>(std::istream& in, DiffClass<T>& v) {
  in >> v.value_;
  int last = in.get();
  if (last != EOF && last != '\t' && last != '\r')
    throw std::invalid_argument("expected tab or EOF");
  return in;
}

std::istream& operator>>(std::istream& in, CategoryClass& v) {
  in >> v.value_;
  int last = in.get();
  if (last != EOF && last != '\t' && last != '\r')
    throw std::invalid_argument("expected tab or EOF");
  return in;
}

typedef DiffClass<int> DaysSinceEpoch;
typedef CategoryClass EventCode; // TODO needs custom distance
typedef CategoryClass QuadClass;
typedef DiffClass<double> GoldsteinScale;
typedef CategoryClass Actor1Name;
typedef CategoryClass Actor2Name;
typedef DiffClass<bool> IsRootEvent;
typedef DiffClass<int> NumMentions;
typedef DiffClass<int> NumSources;
typedef DiffClass<int> NumArticles;
typedef DiffClass<double> AvgTone;
typedef CategoryClass Actor1Geo;
typedef CategoryClass Actor2Geo;
typedef CategoryClass BaseURL;

typedef std::tuple<DaysSinceEpoch, EventCode, QuadClass,
                   GoldsteinScale, Actor1Name, Actor2Name,
                   IsRootEvent, NumMentions, NumSources, NumArticles,
                   AvgTone, Actor1Geo, Actor2Geo, BaseURL> GDELTMini;
