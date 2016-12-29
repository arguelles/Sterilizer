#include <iostream>
#include <iomanip>

int main(){
	char c;
	std::cin >> std::noskipws;
	std::cout << std::uppercase << std::hex << std::setfill('0');
	while(std::cin >> c && !std::cin.eof()){
		if(c==' ')
			std::cout << '+';
		else if(c<'*' || c=='+' || c==',' || c=='/' || (c>=':' && c<='@') || (c>='[' && c<='^') || c=='`' || c>'z')
			std::cout << '%' << std::setw(2) << (int)c;
		else
			std::cout << c;
	}
}