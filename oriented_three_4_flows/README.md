#Oriented version of three 3-flows conjecture

##How to run the code
To compile the code just print "make" in your console (you need g++ compiler).

To run this code you need to specify the number of vertices in the snarks and provide the file with snarks in multicode format.
(you can read instructions in README.md file in three_four_flows folder)

##Main result here (September 18, 2016)
Looks like oriented 23k-flows don't exist for any snark (but I don't have any proof of intuition about why this is so).

On the bright side, seems like oriented 334-flows and oriented 244-flows do always exist.

Exceptions for oriented 333-flows:
* 10 vertices: none
* 18 vertices: none
* 20 vertices: g1, g6
* 22 vertices: g3, g7, g11, g12, g14, g19, g20
* 24 vertices: g1, g3, g16, g17, g21, g24, g25, g27, g28, g30, g31, g35, g36, g37

Oriented abcd- flows:
* Upto and including 22 vertex snarks, there do exist oriented 2233-flows.
