# remember to run chmod 755 ./install.sh first time
# then run ./install.sh
echo ' ______________________________________________________ '
echo '|                                                      |'
echo '|                                                      |'
echo '|               installation begins  !                 |'
echo '|                                                      |'
echo '|______________________________________________________|'                                                
echo
echo 'Removing old buid'
rm -f ./BriXstandalone
echo 'Starting build'
g++ -o BriXstandalone src/main.cpp -std=c++11 || { 
echo "----- Build failed -----"
exit
}
echo 'Build complete'
echo 'Running BriX:'
echo '______________________________________________________'
echo 
./BriXstandalone 