#include "JsonLite.H"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace Foam
{
namespace json
{
namespace
{

void appendEscaped(std::string& out, const std::string& text)
{
    out.push_back('"');
    for (unsigned char c : text)
    {
        switch (c)
        {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b";  break;
            case '\f': out += "\\f";  break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                if (c < 0x20)
                {
                    char buffer[8];
                    std::snprintf(buffer, sizeof(buffer), "\\u%04x",
                                  static_cast<unsigned>(c));
                    out += buffer;
                }
                else
                {
                    // UTF-8 bytes pass through: the one invariant is that the
                    // caller's string is already valid UTF-8 (argv and file
                    // paths on Linux are).
                    out.push_back(static_cast<char>(c));
                }
                break;
        }
    }
    out.push_back('"');
}

void appendDouble(std::string& out, double value)
{
    // JSON has no NaN/Inf. Emitting `null` is the honest choice: a reader must
    // not be able to mistake "not computed" for a measurement.
    if (!std::isfinite(value)) { out += "null"; return; }
    char buffer[40];
    std::snprintf(buffer, sizeof(buffer), "%.17g", value);
    out += buffer;
}

void appendIndent(std::string& out, int indent)
{
    out.append(static_cast<std::size_t>(indent) * 2, ' ');
}

void dumpValue(const Value& v, std::string& out, bool pretty, int depth)
{
    switch (v.type())
    {
        case Value::Type::Null:   out += "null"; return;
        case Value::Type::Bool:   out += v.asBool() ? "true" : "false"; return;
        case Value::Type::Int:
        {
            char buffer[32];
            std::snprintf(buffer, sizeof(buffer), "%lld", v.asInt());
            out += buffer;
            return;
        }
        case Value::Type::Double: appendDouble(out, v.asDouble()); return;
        case Value::Type::String: appendEscaped(out, v.asString()); return;
        case Value::Type::Array:
        {
            if (v.elements().empty()) { out += "[]"; return; }
            out.push_back('[');
            for (std::size_t i = 0; i < v.elements().size(); ++i)
            {
                if (i) out.push_back(',');
                if (pretty) { out.push_back('\n'); appendIndent(out, depth + 1); }
                dumpValue(v.elements()[i], out, pretty, depth + 1);
            }
            if (pretty) { out.push_back('\n'); appendIndent(out, depth); }
            out.push_back(']');
            return;
        }
        case Value::Type::Object:
        {
            if (v.members().empty()) { out += "{}"; return; }
            out.push_back('{');
            for (std::size_t i = 0; i < v.members().size(); ++i)
            {
                if (i) out.push_back(',');
                if (pretty) { out.push_back('\n'); appendIndent(out, depth + 1); }
                appendEscaped(out, v.members()[i].first);
                out.push_back(':');
                if (pretty) out.push_back(' ');
                dumpValue(v.members()[i].second, out, pretty, depth + 1);
            }
            if (pretty) { out.push_back('\n'); appendIndent(out, depth); }
            out.push_back('}');
            return;
        }
    }
}

struct Parser
{
    const char* begin;
    const char* cursor;
    const char* end;
    std::string error;

    void skipSpace()
    {
        while (cursor < end)
        {
            const char c = *cursor;
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') ++cursor;
            else break;
        }
    }

    bool fail(const char* message)
    {
        if (error.empty())
        {
            char buffer[128];
            std::snprintf(buffer, sizeof(buffer), "%s at byte %lld", message,
                          static_cast<long long>(cursor - begin));
            error = buffer;
        }
        return false;
    }

    bool parseValue(Value& out, int depth)
    {
        if (depth > 64) return fail("nesting too deep");
        skipSpace();
        if (cursor >= end) return fail("unexpected end of input");
        const char c = *cursor;
        if (c == '{') return parseObject(out, depth);
        if (c == '[') return parseArray(out, depth);
        if (c == '"') { std::string s; if (!parseString(s)) return false; out = Value(std::move(s)); return true; }
        if (c == 't') { if (!literal("true")) return false; out = Value(true); return true; }
        if (c == 'f') { if (!literal("false")) return false; out = Value(false); return true; }
        if (c == 'n') { if (!literal("null")) return false; out = Value(); return true; }
        return parseNumber(out);
    }

    bool literal(const char* text)
    {
        const std::size_t n = std::strlen(text);
        if (static_cast<std::size_t>(end - cursor) < n) return fail("truncated literal");
        if (std::memcmp(cursor, text, n) != 0) return fail("invalid literal");
        cursor += n;
        return true;
    }

    bool parseNumber(Value& out)
    {
        const char* start = cursor;
        if (cursor < end && (*cursor == '-' || *cursor == '+')) ++cursor;
        bool seenDigit = false, seenDot = false, seenExp = false;
        while (cursor < end)
        {
            const char c = *cursor;
            if (c >= '0' && c <= '9') { seenDigit = true; ++cursor; continue; }
            if (c == '.' && !seenDot && !seenExp) { seenDot = true; ++cursor; continue; }
            if ((c == 'e' || c == 'E') && seenDigit && !seenExp)
            {
                seenExp = true; ++cursor;
                if (cursor < end && (*cursor == '-' || *cursor == '+')) ++cursor;
                continue;
            }
            break;
        }
        if (!seenDigit) return fail("invalid number");
        const std::string token(start, cursor);
        if (!seenDot && !seenExp)
        {
            out = Value(static_cast<long long>(std::strtoll(token.c_str(), nullptr, 10)));
            return true;
        }
        out = Value(std::strtod(token.c_str(), nullptr));
        return true;
    }

    bool parseString(std::string& out)
    {
        if (cursor >= end || *cursor != '"') return fail("expected string");
        ++cursor;
        out.clear();
        while (cursor < end)
        {
            const unsigned char c = static_cast<unsigned char>(*cursor++);
            if (c == '"') return true;
            if (c != '\\') { out.push_back(static_cast<char>(c)); continue; }
            if (cursor >= end) return fail("truncated escape");
            const char e = *cursor++;
            switch (e)
            {
                case '"':  out.push_back('"');  break;
                case '\\': out.push_back('\\'); break;
                case '/':  out.push_back('/');  break;
                case 'b':  out.push_back('\b'); break;
                case 'f':  out.push_back('\f'); break;
                case 'n':  out.push_back('\n'); break;
                case 'r':  out.push_back('\r'); break;
                case 't':  out.push_back('\t'); break;
                case 'u':
                {
                    unsigned code = 0;
                    for (int i = 0; i < 4; ++i)
                    {
                        if (cursor >= end) return fail("truncated \\u escape");
                        const char h = *cursor++;
                        code <<= 4;
                        if (h >= '0' && h <= '9') code |= static_cast<unsigned>(h - '0');
                        else if (h >= 'a' && h <= 'f') code |= static_cast<unsigned>(h - 'a' + 10);
                        else if (h >= 'A' && h <= 'F') code |= static_cast<unsigned>(h - 'A' + 10);
                        else return fail("invalid \\u escape");
                    }
                    // Surrogate pairs are not produced by this writer; encode
                    // the BMP code point and leave astral planes to the reader.
                    if (code < 0x80) out.push_back(static_cast<char>(code));
                    else if (code < 0x800)
                    {
                        out.push_back(static_cast<char>(0xC0 | (code >> 6)));
                        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    }
                    else
                    {
                        out.push_back(static_cast<char>(0xE0 | (code >> 12)));
                        out.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
                        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    }
                    break;
                }
                default: return fail("invalid escape");
            }
        }
        return fail("unterminated string");
    }

    bool parseArray(Value& out, int depth)
    {
        ++cursor; // '['
        out = Value::array();
        skipSpace();
        if (cursor < end && *cursor == ']') { ++cursor; return true; }
        for (;;)
        {
            Value element;
            if (!parseValue(element, depth + 1)) return false;
            out.push(std::move(element));
            skipSpace();
            if (cursor >= end) return fail("unterminated array");
            if (*cursor == ',') { ++cursor; continue; }
            if (*cursor == ']') { ++cursor; return true; }
            return fail("expected ',' or ']'");
        }
    }

    bool parseObject(Value& out, int depth)
    {
        ++cursor; // '{'
        out = Value::object();
        skipSpace();
        if (cursor < end && *cursor == '}') { ++cursor; return true; }
        for (;;)
        {
            skipSpace();
            std::string key;
            if (!parseString(key)) return false;
            skipSpace();
            if (cursor >= end || *cursor != ':') return fail("expected ':'");
            ++cursor;
            Value element;
            if (!parseValue(element, depth + 1)) return false;
            out.set(key, std::move(element));
            skipSpace();
            if (cursor >= end) return fail("unterminated object");
            if (*cursor == ',') { ++cursor; continue; }
            if (*cursor == '}') { ++cursor; return true; }
            return fail("expected ',' or '}'");
        }
    }
};

} // End anonymous namespace

Value& Value::set(const std::string& key, Value v)
{
    if (type_ != Type::Object)
    {
        type_ = Type::Object;
        array_.clear();
    }
    for (auto& member : object_)
        if (member.first == key) { member.second = std::move(v); return member.second; }
    object_.emplace_back(key, std::move(v));
    return object_.back().second;
}

const Value* Value::find(const std::string& key) const
{
    if (type_ != Type::Object) return nullptr;
    for (const auto& member : object_)
        if (member.first == key) return &member.second;
    return nullptr;
}

namespace
{
const Value& nullValue()
{
    static const Value v;
    return v;
}
} // End anonymous namespace

const Value& Value::operator[](const std::string& key) const
{
    const Value* found = find(key);
    return found ? *found : nullValue();
}

Value& Value::operator[](const std::string& key)
{
    if (type_ != Type::Object)
    {
        type_ = Type::Object;
        array_.clear();
        string_.clear();
    }
    for (auto& member : object_)
        if (member.first == key) return member.second;
    object_.emplace_back(key, Value());
    return object_.back().second;
}

Value& Value::push(Value v)
{
    if (type_ != Type::Array)
    {
        type_ = Type::Array;
        object_.clear();
    }
    array_.push_back(std::move(v));
    return array_.back();
}

std::size_t Value::size() const
{
    if (type_ == Type::Array) return array_.size();
    if (type_ == Type::Object) return object_.size();
    return 0;
}

const Value& Value::at(std::size_t i) const
{
    if (type_ != Type::Array || i >= array_.size()) return nullValue();
    return array_[i];
}

std::string Value::dump(bool pretty, int indent) const
{
    std::string out;
    out.reserve(256);
    dumpValue(*this, out, pretty, indent);
    return out;
}

bool parse(const std::string& text, Value& out, std::string& error)
{
    Parser parser{text.data(), text.data(), text.data() + text.size(), std::string()};
    if (!parser.parseValue(out, 0)) { error = parser.error; return false; }
    parser.skipSpace();
    if (parser.cursor != parser.end)
    {
        error = "trailing characters after top-level value";
        return false;
    }
    return true;
}

} // End namespace json
} // End namespace Foam
