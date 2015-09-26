#include "dok.h"

int main(int argc, char const *argv[])
{
	DictionaryOfKeys dic = DictionaryOfKeys(3,4);
	cout << dic(1,2) << endl;
	dic(1,2) = 3.14;
	cout << dic(1,2) << endl;
	dic(1,3) = 4.14;
	cout << dic(1,2) << endl;
	cout << dic(1,3) << endl;
	dic(1,2)= 1.09;
	dic(0,0)= 91.7;
	cout << "=====" << endl;
	cout << dic(1,2) << endl;
	cout << dic(1,3) << endl;
	cout << dic(0,0) << endl;
	dic.Show();
	dic.ShowOctave();
	return 0;
}