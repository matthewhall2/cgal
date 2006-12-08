#!/bin/bash

# Intensive test: test all surface parameterization methods with 1 model in data folder
# Usage: test_model.sh source-file-root
# Example: test_model.sh rotor

echo ""
echo "                    ************************"
echo ""

./test.sh barycentric square taucs eps "$1"
echo "                                -"
<<<<<<< .mine
./test.sh floater circle taucs obj "$1"
=======
./test.sh barycentric circle taucs obj "$1"
echo "                                -"
./test.sh floater circle opennl obj "$1"
>>>>>>> .r35474
echo "                                -"
<<<<<<< .mine
./test.sh conformal square taucs obj "$1"
=======
./test.sh floater square taucs eps "$1"
echo "                                -"
./test.sh conformal circle taucs obj "$1"
>>>>>>> .r35474
echo "                                -"
<<<<<<< .mine
./test.sh authalic circle taucs obj "$1"
=======
./test.sh conformal square opennl eps "$1"
echo "                                -"
./test.sh authalic square taucs obj "$1"
>>>>>>> .r35474
echo "                                -"
<<<<<<< .mine
./test.sh lscm 2pts taucs obj "$1"
=======
./test.sh authalic circle opennl eps "$1"
echo "                                -"
./test.sh lscm 2pts opennl obj "$1"
echo "                                -"
./test.sh lscm 2pts taucs eps "$1"
>>>>>>> .r35474

echo ""
echo "                    ************************"
echo ""

