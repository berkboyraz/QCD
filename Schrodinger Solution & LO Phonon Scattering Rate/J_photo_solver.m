function J_photo= J_photo_solver(ni_matrix, Tijp_matrix, Tipj_matrix, e)

i=1;
j=1;
sum_value=0;
while i<6
    while j<6
        sum_value=ni_matrix(i)/Tijp_matrix(i,j)-ni_matrix(j)/Tipj_matrix(j,i);
        j=j+1;
    end
i=i+1;
end


J_photo = e*sum_value;




end