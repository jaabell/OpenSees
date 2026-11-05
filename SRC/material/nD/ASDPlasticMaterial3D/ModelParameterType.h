// ============================================================================
// Model Parameters base template struct
template<typename T, class NAMER>
struct ModelParameterType {

    static constexpr const char* NAME = NAMER::name;
    T value;
    ModelParameterType() : value() {}   // a parameter not given in Begin_Model_Parameters is 0, not garbage
    ModelParameterType(T x) : value(x) {}
    inline const char* getName() const { return NAME; }
};

//Templated stream insertion
template<typename T, class NAMER>
std::ostream& operator<<(std::ostream& os, const ModelParameterType<T, NAMER>& param) {
    os << param.NAME << " = " << param.value;
    return os;
}
