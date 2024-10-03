using System;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.Xml;
using System.Windows.Forms;

public class Noise
{
    //Pseudorandom Noise functions 1D, 2D, and 3D
    public static double random(double x)
    {
        //pseudorandom 1D noise function, intakes 1 double(Only uses integer part), uses bitwise
        //operators for fast change of the number state, till output as vastly different number
        //multiplied by constant to be decimal
    
        long seedX = (long)(x + 1);
    
        //Constants in each operatiom(21, 11, 37, etc.) were specifically calculated and chosen
        //for calculated and chosen percieved randomness
        long state = seedX ^ (seedX << 21);
        state ^= state << 11;
        state ^= state >> 37;
        state ^= state << 5;
    
        state ^= state >> 20;
        state ^= (state << 17);
        state ^= state >> (int)(seedX % 63 + 1);
    
        //Constant shifts output to be between 0 and 1
        double randomVal = ((double)state) * 0.00000000001;
        randomVal = randomVal - ((int)randomVal);
        if (randomVal < 0)
        {
            randomVal *= -1;
        }
        return randomVal;
    }
    public static double random(double x, double z)
    {
        //pseudorandom 2D noise function, intakes 2 doubles(Only uses integer part), uses bitwise
        //operators for fast change of the number state, till output as vastly different number
        //multiplied by constant to be decimal
    
        long seedX = (long)(x + 1);
        long seedZ = (long)(z + 1);
    
        //Constants in each operatiom(21, 19, 11, 37, etc.) were specifically calculated and chosen
        //for maximum percieved randomness
        long state = seedX ^ (seedX << 21) ^ seedZ ^ (seedZ << 19);
        state ^= state << 11;
        state ^= state >> 37;
        state ^= state << 5;
    
        state ^= state >> 20;
        state ^= (state << 17);
        state ^= state >> (int)(seedZ % 63 + 1);
        state ^= state << (int)(seedX % 63 + 1);
    
        //Constant shifts output to be between 0 and 1
        double randomVal = ((double)state) * 0.00000000001;
        randomVal = randomVal - ((int)randomVal);
        if (randomVal < 0)
        {
            randomVal *= -1;
        }
        return randomVal;
    }
    public static double random(double x, double z, double y)
    {
        //pseudorandom 3D noise function, intakes 3 doubles(Only uses integer part), uses bitwise
        //operators for fast change of the number state, till output as vastly different number
        //multiplied by constant to be decimal
    
        long seedX = (long)(x + 1);
        long seedZ = (long)(z + 1);
        long seedY = (long)(y + 1);
    
        //Constants in each operatiom(21, 19, 23, 11, etc.) were specifically calculated and chosen
        //for maximum percieved randomness
        long state = seedX ^ (seedX << 21) ^ seedZ ^ (seedZ << 19) ^ seedY ^ (seedY << 23);
        state ^= state << 11;
        state ^= state >> 37;
        state ^= state << 5;
    
        state ^= state >> 20;
        state ^= (state << 17);
        state ^= state >> (int)(seedZ % 63 + 1);
        state ^= state << (int)(seedX % 63 + 1);
        state ^= state >> (int)(seedY % 63 + 1);
    
        //Constant shifts output to be between 0 and 1
        double randomVal = ((double)state) * 0.00000000001;
        randomVal = randomVal - ((int)randomVal);
        if (randomVal < 0)
        {
            randomVal *= -1;
        }
        return randomVal;
    }

    //Voronoi Noise Functions
    public static double voronoi(double x, double z)
    {
        //Mirrors this noise from the negative planes (Looks bad where mirrored so be careful)
        x = Math.Abs(x);
        z = Math.Abs(z);
    
        //Put coordinates in a cell
        int cellSize = 200;
        double cellX = (int)(x / cellSize);
        double cellZ = (int)(z / cellSize);
    
        //Run through the 3 by 3 cell area of focus points
        //double focusCellX = cellX;
        //double focusCellZ = cellZ;
        double shortestDist = cellSize * 2;
        double shortestDist2 = shortestDist * 2;
    
        for (int i = -1; i <= 1; i++)
        {
            for (int j = -1; j <= 1; j++)
            {
                double currCellX = i + cellX;
                double currCellZ = j + cellZ;
    
                double xDistance = x - (random(currCellX, currCellZ, 0) * cellSize + cellSize * currCellX);
                double zDistance = z - (random(currCellX, currCellZ, 1) * cellSize + cellSize * currCellZ);
                double distance = Math.Sqrt(xDistance * xDistance + zDistance * zDistance);
    
                //Find closest and 2nd closest nodes to current checked position
                if (distance < shortestDist)
                {
                    shortestDist2 = shortestDist;
                    shortestDist = distance;
                }
                else if (distance < shortestDist2)
                {
                    shortestDist2 = distance;
                }
            }
        }
    
        //Return value between 0 and 1 comparing relationship between closest and 2nd closest nodes
        //return random(focusCellX, focusCellZ); //Returns cells as chunks
        return (shortestDist) / (shortestDist2); //Returns cell boundaries with bounds as bright
    }
    public static double voronoiSecond(double x, double z)
    {
        //Mirrors this noise from the negative planes (Looks bad where mirrored so be careful)
        x = Math.Abs(x);
        z = Math.Abs(z);
    
        //Put coordinates in a cell
        int cellSize = 80;
        double shortestDist = 160;
        double cellX = (int)(x / cellSize);
        double cellZ = (int)(z / cellSize);
    
        //Run through the 3 by 3 cell area of focus points
        double shortestDist2 = shortestDist;
    
        for (int i = -1; i <= 1; i++)
        {
            for (int j = -1; j <= 1; j++)
            {
                double currCellX = i + cellX;
                double currCellZ = j + cellZ;
    
                double xDistance = x - (random(currCellX, currCellZ, 0) * cellSize + cellSize * currCellX);
                double zDistance = z - (random(currCellX, currCellZ, 1) * cellSize + cellSize * currCellZ);
                double distance = Math.Sqrt(xDistance * xDistance + zDistance * zDistance);
    
                //Find closest and 2nd closest nodes to current checked position
                if (distance < shortestDist)
                {
                    shortestDist2 = shortestDist;
                    shortestDist = distance;
                }
                else if (distance < shortestDist2)
                {
                    shortestDist2 = distance;
                }
            }
        }
    
        //Return distance to 2nd closest node and transforming it within a range of 0 to 1
        double maxDist = Math.Sqrt(2) * cellSize;
        return shortestDist2 / maxDist; //Returns cell boundaries with bounds as bright
    }

    //Perlin Noise Functions
    public static double perlin(double x, double z)
    {
        //Mirrors this noise from the negative planes (Sharp flip where mirrored so be careful)
        if (x < 0)
        {
            x = -x;
        }
        if (z < 0)
        {
            z = -z;
        }
    
        //Put coordinates in a cell
        int cellSize = 40; //Chosen Sampling Size
        double baseX = (int)(x / cellSize);
        double baseZ = (int)(z / cellSize);
        double dx = (x % cellSize) / cellSize;
        double dz = (z % cellSize) / cellSize;
    
        //Can check Corners of cells
        /*if(dx < 0.03 && dz < 0.03)
        {
            return 1.0;
        }*/
    
        //Calculate corner vectors
        double[] corners = new double[8]; //Using 1D array for efficiency
        for (int n = 0; n < 4; n++)
        {
            int i = n / 2;
            int j = n % 2;
            double cornerX = baseX + i;
            double cornerZ = baseZ + j;
            double randomVal = random(cornerX, cornerZ);
            randomVal = randomVal * 2 - 1;
            double flip = Math.Sqrt(1 - (randomVal * randomVal)); //Make unit vector
            flip *= ((random(cornerZ, cornerX) > 0.5) ? -1 : 1); //"Randomly" flip direction of unit vector
            corners[n] = randomVal; //Set Vector
            corners[n + 4] = flip;
        }
    
        //Calculate Corner values, calculates with dot product formula
        double a = corners[0] * dx + corners[4] * dz;
        double b = corners[2] * (dx - 1) + corners[6] * dz;
        double c = corners[3] * (dx - 1) + corners[7] * (dz - 1);
        double d = corners[1] * dx + corners[5] * (dz - 1);
    
        //Interpolate between corner values
        double val = interpolation(a, b, c, d, dx, dz);
    
        return (val + 0.71) * 0.7; //Set 0<val<1 (Rounded, true original range +/-sqrt(1/2)=0.7071)
    }
    private static double interpolation(double a, double b, double c, double d, double tx, double tz)
    {
        //Takes inputs and returns one value gotten from weighting a, b, c, and d according to tx and tz
        tx = ((6 * tx - 15) * tx + 10) * tx * tx * tx;
        tz = ((6 * tz - 15) * tz + 10) * tz * tz * tz;//Polynomial interpolation, always t = [0,1), only adjust steepness
        //6t^5 - 15t^4 + 10t^3 = y; 30x^4 - 60x^3 + 30x^2 = dy
        return tx * (-a + b) + tz * (-a + d) + tx * tz * (a - b - d + c) + a;
        //return a*(1 - tx - tz + txz) + b*(tx - txz) + d*(tz - txz) + c*(txz);
        // a - tx a - tz a + txz a + b tx - b txz + d tz - d txz + c txz
        //tx *(-a+b) + tz *(-a+d) + tx * tz *(a-b-d+c) + a
    }
}
