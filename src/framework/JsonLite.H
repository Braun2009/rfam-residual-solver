#ifndef JsonLite_H
#define JsonLite_H

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace Foam
{
namespace json
{

// Minimal JSON document model used by the run telemetry layer.
//
// Why not a third-party parser: the repository links nothing but CoolProp and
// SUNDIALS, and a dependency added for a log format is a dependency added to
// every build. This header implements exactly the subset the telemetry files
// need -- parse, build, serialise -- with insertion-ordered objects so a
// manifest diff between two runs is a textual diff.
//
// Object member order is preserved deliberately: `effective.cfl` must appear
// before `effective.mechanism` in every manifest, forever, or reviewers start
// diffing by eye anyway.
class Value
{
public:
    enum class Type { Null, Bool, Int, Double, String, Array, Object };

    Value() = default;
    Value(std::nullptr_t) {}
    Value(bool v) : type_(Type::Bool), bool_(v) {}
    Value(int v) : type_(Type::Int), int_(v) {}
    Value(long v) : type_(Type::Int), int_(v) {}
    Value(long long v) : type_(Type::Int), int_(v) {}
    Value(unsigned v) : type_(Type::Int), int_(static_cast<long long>(v)) {}
    Value(unsigned long v) : type_(Type::Int), int_(static_cast<long long>(v)) {}
    Value(unsigned long long v) : type_(Type::Int), int_(static_cast<long long>(v)) {}
    Value(double v) : type_(Type::Double), double_(v) {}
    Value(const char* v) : type_(Type::String), string_(v ? v : "") {}
    Value(std::string v) : type_(Type::String), string_(std::move(v)) {}

    static Value array() { Value v; v.type_ = Type::Array; return v; }
    static Value object() { Value v; v.type_ = Type::Object; return v; }

    Type type() const { return type_; }
    bool isNull() const { return type_ == Type::Null; }
    bool isObject() const { return type_ == Type::Object; }
    bool isArray() const { return type_ == Type::Array; }
    bool isString() const { return type_ == Type::String; }
    bool isNumber() const { return type_ == Type::Int || type_ == Type::Double; }
    bool isBool() const { return type_ == Type::Bool; }

    bool asBool(bool fallback = false) const
    {
        if (type_ == Type::Bool) return bool_;
        if (type_ == Type::Int) return int_ != 0;
        return fallback;
    }
    long long asInt(long long fallback = 0) const
    {
        if (type_ == Type::Int) return int_;
        if (type_ == Type::Double) return static_cast<long long>(double_);
        if (type_ == Type::Bool) return bool_ ? 1 : 0;
        return fallback;
    }
    double asDouble(double fallback = 0.0) const
    {
        if (type_ == Type::Double) return double_;
        if (type_ == Type::Int) return static_cast<double>(int_);
        return fallback;
    }
    const std::string& asString() const { return string_; }
    std::string asString(const char* fallback) const
    {
        return type_ == Type::String ? string_ : std::string(fallback ? fallback : "");
    }

    // Object access. `set` appends on first use and overwrites in place after
    // that, so building a manifest is order-stable regardless of code path.
    Value& set(const std::string& key, Value v);
    const Value* find(const std::string& key) const;
    const Value& operator[](const std::string& key) const;
    Value& operator[](const std::string& key);   // creates a null member

    // Array access.
    Value& push(Value v);
    std::size_t size() const;
    const Value& at(std::size_t i) const;

    const std::vector<std::pair<std::string, Value> >& members() const { return object_; }
    const std::vector<Value>& elements() const { return array_; }

    std::string dump(bool pretty = false, int indent = 0) const;

private:
    Type type_ = Type::Null;
    bool bool_ = false;
    long long int_ = 0;
    double double_ = 0.0;
    std::string string_;
    std::vector<Value> array_;
    std::vector<std::pair<std::string, Value> > object_;
};

bool parse(const std::string& text, Value& out, std::string& error);

} // End namespace json
} // End namespace Foam

#endif
