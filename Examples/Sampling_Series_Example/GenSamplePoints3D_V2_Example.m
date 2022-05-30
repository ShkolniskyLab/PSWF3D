function [ SamplePoints ] = GenSamplePoints3D_V2_Example( L )
% This code generates the sampling points in Q that correspond to the
% parameter L, in Cartesian form. Note that the points are constructed in
% the following way: for any azimuth \phi \in [0,pi/2), the points with
% azimuth \phi,\phi+pi/2,\phi+pi,\phi+3pi/2 appear on after the other. The
% last 2*L-1 points are of the form [0,0,z]

counter = 1;

for jx = 1:L
    
   for jy = 0:L
       
      for jz = (-L):L
          
         norm2 = norm([jx/L,jy/L,jz/L]);
                   
         if (norm2 < 1)
             
            SamplePoints(counter,:) = [jx/L,jy/L,jz/L];
            counter = counter+1;
            
            rot = [0 -1 0; 1 0 0; 0 0 1];

            SamplePoints(counter,:) = rot*transpose([jx/L,jy/L,jz/L]);
            counter = counter+1;
            
            SamplePoints(counter,:) = rot*rot*transpose([jx/L,jy/L,jz/L]);
            counter = counter+1;        
            
            SamplePoints(counter,:) = rot*rot*rot*transpose([jx/L,jy/L,jz/L]);
            counter = counter+1;            
            
         end   
            
      end
       
   end
    
end

for jz = (-L+1):(L-1)
               
     SamplePoints(counter,:) = [0,0,jz/L];
     counter = counter+1;
                   
end

end

