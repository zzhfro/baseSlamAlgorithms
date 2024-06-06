#include <iostream>

// 定义一个枚举类型
enum Day {
    Sunday=1,
    Monday,
    Tuesday,
    Wednesday,
    Thursday,
    Friday,
    Saturday
};

int main() {
    // 声明一个枚举变量
    Day today = Monday;
    std::cout<<today<<std::endl;
    // 使用switch语句处理枚举变量
    switch (today) {
        case Sunday:
            std::cout << "Today is Sunday." << std::endl;
            break;
        case Monday:
            std::cout << "Today is Monday." << std::endl;
            break;
        case Tuesday:
            std::cout << "Today is Tuesday." << std::endl;
            break;
        case Wednesday:
            std::cout << "Today is Wednesday." << std::endl;
            break;
        case Thursday:
            std::cout << "Today is Thursday." << std::endl;
            break;
        case Friday:
            std::cout << "Today is Friday." << std::endl;
            break;
        case Saturday:
            std::cout << "Today is Saturday." << std::endl;
            break;
        default:
            std::cerr << "Invalid day." << std::endl;
    }

    return 0;
}
