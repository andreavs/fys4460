function [B , C ] = MK_EQSYSTEM (A , X , Y )
% Total no of internal lattice sites
sites = X *( Y -2);
% Allocate space for the nonzero upper diagonals
main_diag = zeros ( sites , 1);
upper_diag1 = zeros ( sites -1 , 1);
upper_diag2 = zeros ( sites -X , 1);
% Calculates the nonzero upper diagonals
main_diag = A ( X +1: X *( Y -1) , 1) + A ( X +1: X *( Y -1) , 2) + A (1: X *( Y -2) , 2) ...
+ A ( X : X *( Y -1) -1 , 1);
upper_diag1 = A ( X +1: X *( Y -1) -1 , 1);
upper_diag2 = A ( X +1: X *( Y -2) , 2);
main_diag ( find ( main_diag ==0)) = 1;
% Constructing B which is symmetric , lower = upper diagonals .
B = sparse ( sites , sites ); % B * u = t
B = - spdiags ( upper_diag1 , -1 , sites , sites );
B = B + - spdiags ( upper_diag2 , -X , sites , sites );
B = B + B' + spdiags ( main_diag , 0 , sites , sites );
% Constructing C
C = sparse ( sites , 1);
C (1: X ) = A (1: X , 2);
C ( end - X +1: end ) = 0* A (( end -2* X +1: end - X ) , 2);