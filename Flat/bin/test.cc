#include <iostream>
#include <vector>
#include <memory>

using namespace std;


int main (int argc, char ** argv) 
{
    int a = 5;
    shared_ptr<int*> p(make_shared<int*>()); 
    *p = &a;

    shared_ptr<const int*> cp = const_pointer_cast<const int*>(p);

    cout << (**p) << " ";
    a++;
    cout << (**cp) << endl;
}
