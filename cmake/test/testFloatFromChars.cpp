// CMake feature test for floating-point std::from_chars() support

#include <charconv>
#include <string_view>

int main()
{
    const auto s = std::string_view { "2.71828" };
    auto e = 0.0;
    std::from_chars(s.data(), s.data() + s.size(), e);
}
