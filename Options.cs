using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Homework6;

public class Options
{

    // function for Box-Muller
    static Tuple<double, double> BoxMuller()
    {
        Random r = new Random();
        double a = r.NextDouble();
        double b = r.NextDouble();

        double express = Math.Sqrt(-2 * Math.Log(a));

        return Tuple.Create(express * Math.Cos(2 * Math.PI * b), express * Math.Sin(2 * Math.PI * b));
    }

    // generate Gaussian random numbers in matrix given number of paths and time steps
    public static double[,] GenRand(int N, int M)
    {
        double[,] RandNums = new double[M, N];

        for (int j = 0; j < M; j++)
        {
            for (int i = 0; i < N; i++)
            {
                Tuple<double, double> ep = BoxMuller();
                RandNums[j, i] = ep.Item1;
            }
        }
        return RandNums;
    }

    public static Dictionary<string, double> get_outputs(Func<double, double, double, double, double, int, int, double[,], List<double>> func,
                                           double S,
                                           double K,
                                           double T,
                                           double r,
                                           double div,
                                           double sig,
                                           int N,
                                           int M,
                                           double[,] RandNums,
                                           bool IsCall
                                        )
    {
        double dS = 0.001 * S;
        double dsig = 0.001 * sig;
        double dT = 0.001 * T;
        double dr = 0.001 * r;

        Dictionary<string, double> myoutputs = new Dictionary<string, double>();

        // option price
        myoutputs.Add("Price", func(S, T, r, div, sig, N, M, RandNums)[0]);

        // standard error
        myoutputs.Add("Standard Error", func(S, T, r, div, sig, N, M, RandNums)[1]);

        // delta
        myoutputs.Add("Delta", (func(S + dS, T, r, div, sig, N, M, RandNums)[0] - func(S - dS, T, r, div, sig, N, M, RandNums)[0]) / (2 * dS));

        // gamma
        myoutputs.Add("Gamma", (func(S + dS, T, r, div, sig, N, M, RandNums)[0] - 2 * func(S, T, r, div, sig, N, M, RandNums)[0] + func(S - dS, T, r, div, sig, N, M, RandNums)[0]) / (Math.Pow(dS, 2)));
        
        // vega / 100
        myoutputs.Add("Vega / 100", ((func(S, T, r, div, sig + dsig, N, M, RandNums)[0] - func(S, T, r, div, sig - dsig, N, M, RandNums)[0]) / (2 * dsig))/100);

        // theta / 365
        myoutputs.Add("Theta / 365", -((func(S, T + dT, r, div, sig, N, M, RandNums)[0] - func(S, T, r, div, sig, N, M, RandNums)[0]) / dT)/365);
        
        // rho / 100
        myoutputs.Add("Rho / 100", ((func(S, T, r + dr, div, sig, N, M, RandNums)[0] - func(S, T, r - dr, div, sig, N, M, RandNums)[0]) / (2 * dr))/100);

        return myoutputs;
    }

    // function for normal CDF
    // taken from https://jamesmccaffrey.wordpress.com/2014/07/16/the-normal-cumulative-density-function-using-c/
    static double cdf(double z)
    {
        double p = 0.3275911;
        double a1 = 0.254829592;
        double a2 = -0.284496736;
        double a3 = 1.421413741;
        double a4 = -1.453152027;
        double a5 = 1.061405429;

        int sign;
        if (z < 0.0)
            sign = -1;
        else
            sign = 1;

        double x = Math.Abs(z) / Math.Sqrt(2.0);
        double t = 1.0 / (1.0 + p * x);
        double erf = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);
        return 0.5 * (1.0 + sign * erf);
    }

    // function to get delta from Black-Scholes
    static double BlackScholesDelta(double d1, bool IsCall)
    {
        return IsCall is true ? cdf(d1) : cdf(d1) - 1;
    }

    // public enum OptionType
    // {
    //     European,
    //     Asian,
    //     Digital,
    //     Lookback,
    //     Range,
    //     Barrier

    // }

    public class UserInputs
    {
        public double Strike {get; set;}

        public bool is_Call {get; set;}

        public double S {get; set;}

        public double T {get; set;}

        public double r {get; set;}

        public double div {get; set;}

        public double sig {get; set;}

        public int N {get; set;}

        public int M {get; set;}

        public bool antithetic {get; set;}

        public bool multithread {get; set;}

        public bool cont_var {get; set;}

    }

    public class UserInputsBarrier : UserInputs
    {
        public double BarrierLevel {get; set;}

        public int knock {get; set;}
    }

    public class UserInputsDigital : UserInputs
    {
        public double Payout {get; set;}
    }

    public class UserOutputs
    {
        public double Price {get; set;}

        public double SE {get; set;}

        public double Delta {get; set;}

        public double Gamma {get; set;}

        public double Vega {get; set;}

        public double Theta {get; set;}

        public double Rho {get; set;}

        public UserOutputs(double Price, double SE, double Delta, double Gamma, double Vega, double Theta, double Rho)
        {
            this.Price = Price;
            this.SE = SE;
            this.Delta = Delta;
            this.Gamma = Gamma;
            this.Vega = Vega;
            this.Theta = Theta;
            this.Rho = Rho;
        }
    }
    public class European
    {
        public double Strike {get; set;}

        public bool IsCall {get; set;}

        public European(double Strike, bool IsCall)
        {
            this.Strike = Strike;
            this.IsCall = IsCall;
        }

        public List<double> GetPrice(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double lnSt = lnS;

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                }

                double ST = Math.Exp(lnSt);
                double OT = Math.Max(option_type * (ST - this.Strike), 0);

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double lnSt = lnS;

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];  // get number for random number matrix
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                        }

                        double ST = Math.Exp(lnSt);
                        double OT = Math.Max(option_type * (ST - this.Strike), 0);

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum() * sum_OT.Sum()/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig * Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0; // for standard error calculation

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)
            {
                double lnSt1 = lnS;
                double lnSt2 = lnS;

                for (int i = 0; i < N; i++)
                {
                    lnSt1 = lnSt1 + nudt + sigsdt * RandNums[j,i];
                    lnSt2 = lnSt2 + nudt + sigsdt * -RandNums[j,i];
                }

                double St1 = Math.Exp(lnSt1);
                double St2 = Math.Exp(lnSt2);
                
                double OT = 0.5 * (Math.Max(option_type * (St1 - this.Strike), 0) + Math.Max(option_type * (St2 - this.Strike), 0));

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT*sum_OT/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig * Math.Sqrt(dt);
            double lnS = Math.Log(S);
            
            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double lnSt1 = lnS;
                        double lnSt2 = lnS;

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            lnSt1 = lnSt1 + nudt + sigsdt * RandNums[d,i];
                            lnSt2 = lnSt2 + nudt + sigsdt * -RandNums[d,i];
                        }

                        double St1 = Math.Exp(lnSt1);
                        double St2 = Math.Exp(lnSt2);

                        double OT = 0.5 * (Math.Max(option_type * (St1 - this.Strike), 0) + Math.Max(option_type * (St2 - this.Strike), 0));

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;
            double nudt = (r - div - 0.5 * Math.Pow(sig, 2)) * dt;
            double sigsdt = sig * Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            int beta1 = -1;

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)
            {
                double St = S;
                double cv = 0;

                for (int i = 0; i < N; i++)
                {
                    double t = (i - 1)*dt;
                    double d1 = (Math.Log(St / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double delta = BlackScholesDelta(d1, this.IsCall);
                    double epsilon = RandNums[j,i];
                    double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                    cv = cv + delta * (Stn - St * erddt);
                    St = Stn;
                }
                
                double OT = Math.Max(option_type * (St - this.Strike), 0) + beta1 * cv;
                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r*T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT*sum_OT/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;
            double nudt = (r - div - 0.5 * Math.Pow(sig, 2)) * dt;
            double sigsdt = sig * Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            int beta1 = -1;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {

                        double St = S;
                        double cv = 0;

                        for (int i = 0; i < N; i++)
                        {
                            double t = (i - 1)*dt;
                            double d1 = (Math.Log(St / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double delta = BlackScholesDelta(d1, this.IsCall);
                            double epsilon = RandNums[d,i];
                            double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                            cv = cv + delta * (Stn - St * erddt);
                            St = Stn;
                        }
                
                        double OT = Math.Max(option_type * (St - this.Strike), 0) + beta1 * cv;
                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r*T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;
            double nudt = (r - div - 0.5 * Math.Pow(sig, 2)) * dt;
            double sigsdt = sig * Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            int beta1 = -1;

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)
            {
                double St1 = S;
                double St2 = S;
                double cv1 = 0;
                double cv2 = 0;

                for (int i = 0; i < N; i++)
                {
                    double t = (i - 1) * dt;

                    double d1 = (Math.Log(St1 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double d1_2 = (Math.Log(St2 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    
                    double delta1 = BlackScholesDelta(d1, this.IsCall);
                    double delta2 = BlackScholesDelta(d1_2, this.IsCall);
                    
                    double epsilon = RandNums[j,i];

                    double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                    double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                    cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                    cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                    St1 = Stn1;
                    St2 = Stn2;
                }

                // option terminal
                double OT = 0.5 * (Math.Max(option_type * (St1 - this.Strike), 0) + beta1 * cv1 +
                                    Math.Max(option_type * (St2 - this.Strike), 0) + beta1 * cv2);

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);

            double SD = Math.Sqrt( (sum_OT2 - sum_OT*sum_OT/M) * Math.Exp(-2 * r * T) / (M - 1) );
            double SE = SD / Math.Sqrt(M);

            // return option value and standard error in a list
            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;
            double nudt = (r - div - 0.5 * Math.Pow(sig, 2)) * dt;
            double sigsdt = sig * Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            int beta1 = -1;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double St1 = S;
                        double St2 = S;
                        double cv1 = 0;
                        double cv2 = 0;

                        for (int i = 0; i < N; i++)
                        {
                            double t = (i - 1) * dt;

                            double d1 = (Math.Log(St1 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double d1_2 = (Math.Log(St2 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    
                            double delta1 = BlackScholesDelta(d1, this.IsCall);
                            double delta2 = BlackScholesDelta(d1_2, this.IsCall);
                    
                            double epsilon = RandNums[d,i];

                            double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                            double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                            cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                            cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                            St1 = Stn1;
                            St2 = Stn2;
                        }

                        // option terminal
                        double OT = 0.5 * (Math.Max(option_type * (St1 - this.Strike), 0) + beta1 * cv1 + Math.Max(option_type * (St2 - this.Strike), 0) + beta1 * cv2);

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);

            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M) * Math.Exp(-2 * r * T) / (M - 1) );
            double SE = SD / Math.Sqrt(M);

            // return option value and standard error in a list
            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }
    }

    public class Asian : Option
    {
        public double Strike {get; set;}

        public bool IsCall {get; set;}

        public Asian(double Strike, bool IsCall)
        {
            this.Strike = Strike;
            this.IsCall = IsCall;
        }

        public List<double> GetPrice(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double lnSt = lnS;
                List<double> path_St = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                    path_St.Add(Math.Exp(lnSt));
                }

                double OT = Math.Max(option_type * (path_St.Average() - this.Strike), 0);

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++)
            {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++)
                    {
                        double lnSt = lnS;
                        List<double> path_St = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                            path_St.Add(Math.Exp(lnSt));
                        }

                        double OT = Math.Max(option_type * (path_St.Average() - this.Strike), 0);

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum() * sum_OT.Sum()/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double lnSt = lnS;
                double lnSt2 = lnS;
                List<double> path_St = new List<double>();
                List<double> path_St2 = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                    lnSt2 = lnSt2 + nudt + sigsdt * -epsilon;
                    path_St.Add(Math.Exp(lnSt));
                    path_St2.Add(Math.Exp(lnSt));
                }

                double OT = 0.5 * (Math.Max(option_type * (path_St.Average() - this.Strike), 0) + Math.Max(option_type * (path_St2.Average() - this.Strike), 0));

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double lnSt = lnS;
                        double lnSt2 = lnS;
                        List<double> path_St = new List<double>();
                        List<double> path_St2 = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                            lnSt2 = lnSt2 + nudt + sigsdt * -epsilon;
                            path_St.Add(Math.Exp(lnSt));
                            path_St2.Add(Math.Exp(lnSt));
                        }

                        double OT = 0.5 * (Math.Max(option_type * (path_St.Average() - this.Strike), 0) + Math.Max(option_type * (path_St2.Average() - this.Strike), 0));

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double St = S;
                double cv = 0;
                List<double> path_St = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double t = (i - 1)*dt;
                    double d1 = (Math.Log(St / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double delta = BlackScholesDelta(d1, this.IsCall);
                    double epsilon = RandNums[j,i];
                    double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                    cv = cv + delta * (Stn - St * erddt);
                    St = Stn;
                    path_St.Add(St);
                }

                double OT = Math.Max(option_type * (path_St.Average() - this.Strike), 0) + beta1 * cv;

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {

                        double St = S;
                        double cv = 0;
                        List<double> path_St = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double t = (i - 1)*dt;
                            double d1 = (Math.Log(St / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double delta = BlackScholesDelta(d1, this.IsCall);
                            double epsilon = RandNums[d,i];
                            double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                            cv = cv + delta * (Stn - St * erddt);
                            St = Stn;
                            path_St.Add(St);
                        }

                        double OT = Math.Max(option_type * (path_St.Average() - this.Strike), 0) + beta1 * cv;

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r*T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double St1 = S;
                double St2 = S;
                double cv1 = 0;
                double cv2 = 0;
                List<double> path_St1 = new List<double>();
                List<double> path_St2 = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double t = (i - 1) * dt;

                    double d1 = (Math.Log(St1 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double d1_2 = (Math.Log(St2 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    
                    double delta1 = BlackScholesDelta(d1, this.IsCall);
                    double delta2 = BlackScholesDelta(d1_2, this.IsCall);
                    
                    double epsilon = RandNums[j,i];

                    double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                    double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                    cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                    cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                    St1 = Stn1;
                    St2 = Stn2;

                    path_St1.Add(St1);
                    path_St2.Add(St2);
                }

                double OT = 0.5 * (Math.Max(option_type * (path_St1.Average() - this.Strike), 0) + beta1 * cv1 +
                                    Math.Max(option_type * (path_St2.Average() - this.Strike), 0) + beta1 * cv2);

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double St1 = S;
                        double St2 = S;
                        double cv1 = 0;
                        double cv2 = 0;
                        List<double> path_St1 = new List<double>();
                        List<double> path_St2 = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double t = (i - 1) * dt;

                            double d1 = (Math.Log(St1 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double d1_2 = (Math.Log(St2 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            
                            double delta1 = BlackScholesDelta(d1, this.IsCall);
                            double delta2 = BlackScholesDelta(d1_2, this.IsCall);
                            
                            double epsilon = RandNums[d,i];

                            double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                            double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                            cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                            cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                            St1 = Stn1;
                            St2 = Stn2;

                            path_St1.Add(St1);
                            path_St2.Add(St2);
                        }

                        double OT = 0.5 * (Math.Max(option_type * (path_St1.Average() - this.Strike), 0) + beta1 * cv1 +
                                        Math.Max(option_type * (path_St2.Average() - this.Strike), 0) + beta1 * cv2);

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M) * Math.Exp(-2 * r * T) / (M - 1) );
            double SE = SD / Math.Sqrt(M);

            // return option value and standard error in a list
            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }
    }

    public class Digital : Option
    {
        public double Strike {get; set;}

        public bool IsCall {get; set;}

        public double Payout {get; set;}

        public Digital(double Strike, bool IsCall, double Payout)
        {
            this.Strike = Strike;
            this.IsCall = IsCall;
            this.Payout = Payout;
        }

        public List<double> GetPrice(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            Func<double, double> get_OT = null;

            if (IsCall is true)
                get_OT = x => Convert.ToInt32(x > this.Strike) * this.Payout;
            else
                get_OT = x => Convert.ToInt32(x < this.Strike) * this.Payout;
            
            for (int j = 0; j < M; j++)  // for each simulation
            {
                double lnSt = lnS;

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                }

                double ST = Math.Exp(lnSt);

                double OT = get_OT(ST);        

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);

            return MyOutput;
        }

        public List<double> GetPrice_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            Func<double, double> get_OT = null;

            if (IsCall is true)
                get_OT = x => Convert.ToInt32(x > this.Strike) * this.Payout;
            else
                get_OT = x => Convert.ToInt32(x < this.Strike) * this.Payout;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++)
            {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double lnSt = lnS;

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                        }

                        double ST = Math.Exp(lnSt);

                        double OT = get_OT(ST);        

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum() * sum_OT.Sum()/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            Func<double, double> get_OT = null;

            if (IsCall is true)
                get_OT = x => Convert.ToInt32(x > this.Strike) * this.Payout;
            else
                get_OT = x => Convert.ToInt32(x < this.Strike) * this.Payout;
            
            for (int j = 0; j < M; j++)  // for each simulation
            {
                double lnSt = lnS;
                double lnSt2 = lnS;

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                    lnSt2 = lnSt2 + nudt + sigsdt * -epsilon;
                }

                double ST = Math.Exp(lnSt);
                double ST2 = Math.Exp(lnSt2);

                double OT = 0.5 * (get_OT(ST) + get_OT(ST2));

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);

            return MyOutput;
        }

        public List<double> GetPrice_Anti_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            Func<double, double> get_OT = null;

            if (IsCall is true)
                get_OT = x => Convert.ToInt32(x > this.Strike) * this.Payout;
            else
                get_OT = x => Convert.ToInt32(x < this.Strike) * this.Payout;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];
            
            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double lnSt = lnS;
                        double lnSt2 = lnS;

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                            lnSt2 = lnSt2 + nudt + sigsdt * -epsilon;
                        }

                        double ST = Math.Exp(lnSt);
                        double ST2 = Math.Exp(lnSt2);

                        double OT = 0.5 * (get_OT(ST) + get_OT(ST2));

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            double sum_OT = 0;
            double sum_OT2 = 0;

            Func<double, double> get_OT = null;

            if (IsCall is true)
                get_OT = x => Convert.ToInt32(x > this.Strike) * this.Payout;
            else
                get_OT = x => Convert.ToInt32(x < this.Strike) * this.Payout;
            
            for (int j = 0; j < M; j++)  // for each simulation
            {
                double St = S;
                double cv = 0;

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double t = (i - 1)*dt;
                    double d1 = (Math.Log(St / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double delta = BlackScholesDelta(d1, this.IsCall);
                    double epsilon = RandNums[j,i];
                    double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                    cv = cv + delta * (Stn - St * erddt);
                    St = Stn;
                }

                double OT = get_OT(St) + beta1 * cv;     

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);

            return MyOutput;
        }

        public List<double> GetPrice_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            Func<double, double> get_OT = null;

            if (IsCall is true)
                get_OT = x => Convert.ToInt32(x > this.Strike) * this.Payout;
            else
                get_OT = x => Convert.ToInt32(x < this.Strike) * this.Payout;
            
            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++)
                    {
                        double St = S;
                        double cv = 0;

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double t = (i - 1)*dt;
                            double d1 = (Math.Log(St / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double delta = BlackScholesDelta(d1, this.IsCall);
                            double epsilon = RandNums[d,i];
                            double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                            cv = cv + delta * (Stn - St * erddt);
                            St = Stn;
                        }

                        double OT = get_OT(St) + beta1 * cv;     

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r*T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            double sum_OT = 0;
            double sum_OT2 = 0;

            Func<double, double> get_OT = null;

            if (IsCall is true)
                get_OT = x => Convert.ToInt32(x > this.Strike) * this.Payout;
            else
                get_OT = x => Convert.ToInt32(x < this.Strike) * this.Payout;
            
            for (int j = 0; j < M; j++)  // for each simulation
            {
                double St1 = S;
                double St2 = S;
                double cv1 = 0;
                double cv2 = 0;

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double t = (i - 1) * dt;

                    double d1 = (Math.Log(St1 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double d1_2 = (Math.Log(St2 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    
                    double delta1 = BlackScholesDelta(d1, this.IsCall);
                    double delta2 = BlackScholesDelta(d1_2, this.IsCall);
                    
                    double epsilon = RandNums[j,i];

                    double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                    double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                    cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                    cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                    St1 = Stn1;
                    St2 = Stn2;
                }

                double OT = 0.5 * (get_OT(St1) + beta1 * cv1 +
                                    get_OT(St2) + beta1 * cv2);

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);

            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            Func<double, double> get_OT = null;

            if (IsCall is true)
                get_OT = x => Convert.ToInt32(x > this.Strike) * this.Payout;
            else
                get_OT = x => Convert.ToInt32(x < this.Strike) * this.Payout;
            
            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double St1 = S;
                        double St2 = S;
                        double cv1 = 0;
                        double cv2 = 0;

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double t = (i - 1) * dt;

                            double d1 = (Math.Log(St1 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double d1_2 = (Math.Log(St2 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            
                            double delta1 = BlackScholesDelta(d1, this.IsCall);
                            double delta2 = BlackScholesDelta(d1_2, this.IsCall);
                            
                            double epsilon = RandNums[d,i];

                            double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                            double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                            cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                            cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                            St1 = Stn1;
                            St2 = Stn2;
                        }

                        double OT = 0.5 * (get_OT(St1) + beta1 * cv1 +
                                        get_OT(St2) + beta1 * cv2);

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);

            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M) * Math.Exp(-2 * r * T) / (M - 1) );
            double SE = SD / Math.Sqrt(M);

            // return option value and standard error in a list
            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }
    }

    public class Lookback : Option
    {
        public double Strike {get; set;}

        public bool IsCall {get; set;}

        public Lookback(double Strike, bool IsCall)
        {
            this.Strike = Strike;
            this.IsCall = IsCall;
        }

        public List<double> GetPrice(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            // to be used in the payoff function; 1 if Call; -1 if Put
            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double lnSt = lnS;
                List<double> path_St = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                    path_St.Add(Math.Exp(lnSt));
                }

                double OT = Math.Max(option_type * (path_St.Max() - this.Strike), 0);

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            // to be used in the payoff function; 1 if Call; -1 if Put
            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++)
            {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double lnSt = lnS;
                        List<double> path_St = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                            path_St.Add(Math.Exp(lnSt));
                        }

                        double OT = Math.Max(option_type * (path_St.Max() - this.Strike), 0);

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum() * sum_OT.Sum()/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            // to be used in the payoff function; 1 if Call; -1 if Put
            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double lnSt = lnS;
                double lnSt2 = lnS;
                List<double> path_St = new List<double>();
                List<double> path_St2 = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                    lnSt2 = lnSt2 + nudt + sigsdt * -epsilon;
                    path_St.Add(Math.Exp(lnSt));
                    path_St2.Add(Math.Exp(lnSt2));
                }

                double OT = 0.5 * (Math.Max(option_type * (path_St.Max() - this.Strike), 0) + Math.Max(option_type * (path_St2.Max() - this.Strike), 0));

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            // to be used in the payoff function; 1 if Call; -1 if Put
            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double lnSt = lnS;
                        double lnSt2 = lnS;
                        List<double> path_St = new List<double>();
                        List<double> path_St2 = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                            lnSt2 = lnSt2 + nudt + sigsdt * -epsilon;
                            path_St.Add(Math.Exp(lnSt));
                            path_St2.Add(Math.Exp(lnSt2));
                        }

                        double OT = 0.5 * (Math.Max(option_type * (path_St.Max() - this.Strike), 0) + Math.Max(option_type * (path_St2.Max() - this.Strike), 0));

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            double sum_OT = 0;
            double sum_OT2 = 0;

            // to be used in the payoff function; 1 if Call; -1 if Put
            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double St = S;
                double cv = 0;
                List<double> path_St = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double t = (i - 1)*dt;
                    double d1 = (Math.Log(St / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double delta = BlackScholesDelta(d1, this.IsCall);
                    double epsilon = RandNums[j,i];
                    double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                    cv = cv + delta * (Stn - St * erddt);
                    St = Stn;
                    path_St.Add(St);
                }

                double OT = Math.Max(option_type * (path_St.Max() - this.Strike), 0) + beta1 * cv;

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            // to be used in the payoff function; 1 if Call; -1 if Put
            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++)
                    {
                        double St = S;
                        double cv = 0;
                        List<double> path_St = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double t = (i - 1)*dt;
                            double d1 = (Math.Log(St / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double delta = BlackScholesDelta(d1, this.IsCall);
                            double epsilon = RandNums[d,i];
                            double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                            cv = cv + delta * (Stn - St * erddt);
                            St = Stn;
                            path_St.Add(St);
                        }

                        double OT = Math.Max(option_type * (path_St.Max() - this.Strike), 0) + beta1 * cv;

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r*T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            double sum_OT = 0;
            double sum_OT2 = 0;

            // to be used in the payoff function; 1 if Call; -1 if Put
            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double St1 = S;
                double St2 = S;
                double cv1 = 0;
                double cv2 = 0;
                List<double> path_St1 = new List<double>();
                List<double> path_St2 = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double t = (i - 1) * dt;

                    double d1 = (Math.Log(St1 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double d1_2 = (Math.Log(St2 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    
                    double delta1 = BlackScholesDelta(d1, this.IsCall);
                    double delta2 = BlackScholesDelta(d1_2, this.IsCall);
                    
                    double epsilon = RandNums[j,i];

                    double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                    double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                    cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                    cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                    St1 = Stn1;
                    St2 = Stn2;

                    path_St1.Add(St1);
                    path_St2.Add(St2);
                }

                double OT = 0.5 * (Math.Max(option_type * (path_St1.Max() - this.Strike), 0) + beta1 * cv1 +
                                    Math.Max(option_type * (path_St2.Max() - this.Strike), 0) + beta1 * cv2);

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            // to be used in the payoff function; 1 if Call; -1 if Put
            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double St1 = S;
                        double St2 = S;
                        double cv1 = 0;
                        double cv2 = 0;
                        List<double> path_St1 = new List<double>();
                        List<double> path_St2 = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double t = (i - 1) * dt;

                            double d1 = (Math.Log(St1 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double d1_2 = (Math.Log(St2 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            
                            double delta1 = BlackScholesDelta(d1, this.IsCall);
                            double delta2 = BlackScholesDelta(d1_2, this.IsCall);
                            
                            double epsilon = RandNums[d,i];

                            double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                            double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                            cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                            cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                            St1 = Stn1;
                            St2 = Stn2;

                            path_St1.Add(St1);
                            path_St2.Add(St2);
                        }

                        double OT = 0.5 * (Math.Max(option_type * (path_St1.Max() - this.Strike), 0) + beta1 * cv1 +
                                        Math.Max(option_type * (path_St2.Max() - this.Strike), 0) + beta1 * cv2);

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);

            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M) * Math.Exp(-2 * r * T) / (M - 1) );
            double SE = SD / Math.Sqrt(M);

            // return option value and standard error in a list
            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }
    }

    public class Range : Option
    {
        public List<double> GetPrice(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double lnSt = lnS;
                List<double> path_St = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                    path_St.Add(Math.Exp(lnSt));
                }

                double OT = path_St.Max() - path_St.Min();

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++)
            {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double lnSt = lnS;
                        List<double> path_St = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                            path_St.Add(Math.Exp(lnSt));
                        }

                        double OT = path_St.Max() - path_St.Min();

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum() * sum_OT.Sum()/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double lnSt = lnS;
                double lnSt2 = lnS;
                List<double> path_St = new List<double>();
                List<double> path_St2 = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                    lnSt2 = lnSt2 + nudt + sigsdt * -epsilon;
                    path_St.Add(Math.Exp(lnSt));
                    path_St2.Add(Math.Exp(lnSt2));
                }

                double OT = 0.5 * ((path_St.Max() - path_St.Min()) + (path_St2.Max() - path_St2.Min()));

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double lnSt = lnS;
                        double lnSt2 = lnS;
                        List<double> path_St = new List<double>();
                        List<double> path_St2 = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                            lnSt2 = lnSt2 + nudt + sigsdt * -epsilon;
                            path_St.Add(Math.Exp(lnSt));
                            path_St2.Add(Math.Exp(lnSt2));
                        }

                        double OT = 0.5 * ((path_St.Max() - path_St.Min()) + (path_St2.Max() - path_St2.Min()));

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            double sum_OT = 0;
            double sum_OT2 = 0;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double St = S;
                double cv = 0;
                List<double> path_St = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double t = (i - 1)*dt;
                    double d1 = (Math.Log(St) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double delta = BlackScholesDelta(d1, true);
                    double epsilon = RandNums[j,i];
                    double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                    cv = cv + delta * (Stn - St * erddt);
                    St = Stn;
                    path_St.Add(St);
                }

                double OT = path_St.Max() - path_St.Min() + beta1 * cv;

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++)
                    {
                        double St = S;
                        double cv = 0;
                        List<double> path_St = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double t = (i - 1)*dt;
                            double d1 = (Math.Log(St) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double delta = BlackScholesDelta(d1, true);
                            double epsilon = RandNums[d,i];
                            double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                            cv = cv + delta * (Stn - St * erddt);
                            St = Stn;
                            path_St.Add(St);
                        }

                        double OT = path_St.Max() - path_St.Min() + beta1 * cv;

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r*T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            double sum_OT = 0;
            double sum_OT2 = 0;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double St1 = S;
                double St2 = S;
                double cv1 = 0;
                double cv2 = 0;
                List<double> path_St1 = new List<double>();
                List<double> path_St2 = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double t = (i - 1) * dt;

                    double d1 = (Math.Log(St1) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double d1_2 = (Math.Log(St2) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    
                    double delta1 = BlackScholesDelta(d1, true);
                    double delta2 = BlackScholesDelta(d1_2, true);
                    
                    double epsilon = RandNums[j,i];

                    double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                    double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                    cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                    cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                    St1 = Stn1;
                    St2 = Stn2;
                    path_St1.Add(St1);
                    path_St2.Add(St2);
                }

                double OT = 0.5 * (path_St1.Max() - path_St1.Min() + beta1 * cv1 +
                                    path_St2.Max() - path_St2.Min() + beta1 * cv2);

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++)
                    {
                        double St1 = S;
                        double St2 = S;
                        double cv1 = 0;
                        double cv2 = 0;
                        List<double> path_St1 = new List<double>();
                        List<double> path_St2 = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double t = (i - 1) * dt;

                            double d1 = (Math.Log(St1) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double d1_2 = (Math.Log(St2) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            
                            double delta1 = BlackScholesDelta(d1, true);
                            double delta2 = BlackScholesDelta(d1_2, true);
                            
                            double epsilon = RandNums[d,i];

                            double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                            double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                            cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                            cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                            St1 = Stn1;
                            St2 = Stn2;
                            path_St1.Add(St1);
                            path_St2.Add(St2);
                        }

                        double OT = 0.5 * (path_St1.Max() - path_St1.Min() + beta1 * cv1 +
                                        path_St2.Max() - path_St2.Min() + beta1 * cv2);

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);

            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M) * Math.Exp(-2 * r * T) / (M - 1) );
            double SE = SD / Math.Sqrt(M);

            // return option value and standard error in a list
            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }
    }

    public enum KnockType
    {
        DownAndOut,
        UpAndOut,
        DownAndIn,
        UpAndIn

    }

    public class Barrier : Option
    {
        public double Strike {get; set;}

        public bool IsCall {get; set;}

        public double BarrierLevel {get; set;}

        public KnockType KnockType {get; set;}

        public Barrier(double Strike, bool IsCall, double BarrierLevel, KnockType KnockType)
        {
            this.Strike = Strike;
            this.IsCall = IsCall;
            this.BarrierLevel = BarrierLevel;
            this.KnockType = KnockType;
        }

        public List<double> GetPrice(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double lnSt = lnS;
                List<double> path_St = new List<double>();

                bool barrier_hit = false;
                int flag = 0;

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                    path_St.Add(Math.Exp(lnSt));
                }

                if (this.KnockType == KnockType.DownAndOut)
                    barrier_hit = path_St.Where(x => x < this.BarrierLevel).Count() > 0;
                else if (this.KnockType == KnockType.UpAndOut)
                    barrier_hit = path_St.Where(x => x > this.BarrierLevel).Count() > 0;
                else if (this.KnockType == KnockType.DownAndIn)
                    barrier_hit = path_St.Where(x => x < this.BarrierLevel).Count() > 0;
                else if (this.KnockType == KnockType.UpAndIn)
                    barrier_hit = path_St.Where(x => x > this.BarrierLevel).Count() > 0;

                if (barrier_hit is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag = 0;
                else if (barrier_hit is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag = 1;
                else if (barrier_hit is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag = 1;
                else if (barrier_hit is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag = 0;

                double OT = Math.Max(flag * option_type * (Math.Exp(lnSt) - this.Strike), 0);

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++)
            {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        double lnSt = lnS;
                        List<double> path_St = new List<double>();

                        bool barrier_hit = false;
                        int flag = 0;

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                            path_St.Add(Math.Exp(lnSt));
                        }

                        if (this.KnockType == KnockType.DownAndOut)
                            barrier_hit = path_St.Where(x => x < this.BarrierLevel).Count() > 0;
                        else if (this.KnockType == KnockType.UpAndOut)
                            barrier_hit = path_St.Where(x => x > this.BarrierLevel).Count() > 0;
                        else if (this.KnockType == KnockType.DownAndIn)
                            barrier_hit = path_St.Where(x => x < this.BarrierLevel).Count() > 0;
                        else if (this.KnockType == KnockType.UpAndIn)
                            barrier_hit = path_St.Where(x => x > this.BarrierLevel).Count() > 0;

                        if (barrier_hit is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag = 0;
                        else if (barrier_hit is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag = 1;
                        else if (barrier_hit is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag = 1;
                        else if (barrier_hit is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag = 0;

                        double OT = Math.Max(flag * option_type * (Math.Exp(lnSt) - this.Strike), 0);

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum() * sum_OT.Sum()/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                bool barrier_hit1 = false;
                bool barrier_hit2 = false;
                int flag1 = 0;
                int flag2 = 0;

                double lnSt = lnS;
                double lnSt2 = lnS;
                List<double> path_St1 = new List<double>();
                List<double> path_St2 = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double epsilon = RandNums[j, i];
                    lnSt = lnSt + nudt + sigsdt * epsilon;
                    lnSt2 = lnSt2 + nudt + sigsdt * -epsilon;
                    path_St1.Add(Math.Exp(lnSt));
                    path_St2.Add(Math.Exp(lnSt2));
                }

                if (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.DownAndIn)
                {
                    barrier_hit1 = path_St1.Where(x => x < this.BarrierLevel).Count() > 0;
                    barrier_hit2 = path_St2.Where(x => x < this.BarrierLevel).Count() > 0;
                }
                else if (this.KnockType == KnockType.UpAndOut | this.KnockType == KnockType.UpAndIn)
                {
                    barrier_hit1 = path_St1.Where(x => x > this.BarrierLevel).Count() > 0;
                    barrier_hit2 = path_St2.Where(x => x > this.BarrierLevel).Count() > 0;
                }

                if (barrier_hit1 is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag1 = 0;
                else if (barrier_hit1 is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag1 = 1;
                else if (barrier_hit1 is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag1 = 1;
                else if (barrier_hit1 is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag1 = 0;

                if (barrier_hit2 is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag2 = 0;
                else if (barrier_hit2 is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag2 = 1;
                else if (barrier_hit2 is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag2 = 1;
                else if (barrier_hit2 is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag2 = 0;

                double OT = 0.5 * (Math.Max(flag1 * option_type * (Math.Exp(lnSt) - this.Strike), 0) + Math.Max(flag2 * option_type * (Math.Exp(lnSt2) - this.Strike), 0));

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double lnS = Math.Log(S);

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        bool barrier_hit1 = false;
                        bool barrier_hit2 = false;
                        int flag1 = 0;
                        int flag2 = 0;

                        double lnSt = lnS;
                        double lnSt2 = lnS;
                        List<double> path_St1 = new List<double>();
                        List<double> path_St2 = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double epsilon = RandNums[d, i];
                            lnSt = lnSt + nudt + sigsdt * epsilon;
                            lnSt2 = lnSt2 + nudt + sigsdt * -epsilon;
                            path_St1.Add(Math.Exp(lnSt));
                            path_St2.Add(Math.Exp(lnSt2));
                        }

                        if (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.DownAndIn)
                        {
                            barrier_hit1 = path_St1.Where(x => x < this.BarrierLevel).Count() > 0;
                            barrier_hit2 = path_St2.Where(x => x < this.BarrierLevel).Count() > 0;
                        }
                        else if (this.KnockType == KnockType.UpAndOut | this.KnockType == KnockType.UpAndIn)
                        {
                            barrier_hit1 = path_St1.Where(x => x > this.BarrierLevel).Count() > 0;
                            barrier_hit2 = path_St2.Where(x => x > this.BarrierLevel).Count() > 0;
                        }

                        if (barrier_hit1 is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag1 = 0;
                        else if (barrier_hit1 is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag1 = 1;
                        else if (barrier_hit1 is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag1 = 1;
                        else if (barrier_hit1 is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag1 = 0;

                        if (barrier_hit2 is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag2 = 0;
                        else if (barrier_hit2 is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag2 = 1;
                        else if (barrier_hit2 is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag2 = 1;
                        else if (barrier_hit2 is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag2 = 0;

                        double OT = 0.5 * (Math.Max(flag1 * option_type * (Math.Exp(lnSt) - this.Strike), 0) + Math.Max(flag2 * option_type * (Math.Exp(lnSt2) - this.Strike), 0));

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                double St = S;
                double cv = 0;
                List<double> path_St = new List<double>();

                bool barrier_hit = false;
                int flag = 0;

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double t = (i - 1)*dt;
                    double d1 = (Math.Log(St / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double delta = BlackScholesDelta(d1, this.IsCall);
                    double epsilon = RandNums[j,i];
                    double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                    cv = cv + delta * (Stn - St * erddt);
                    St = Stn;
                    path_St.Add(St);
                }

                if (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.DownAndIn)
                    barrier_hit = path_St.Where(x => x < this.BarrierLevel).Count() > 0;
                else if (this.KnockType == KnockType.UpAndOut | this.KnockType == KnockType.UpAndIn)
                    barrier_hit = path_St.Where(x => x > this.BarrierLevel).Count() > 0;

                if (barrier_hit is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag = 0;
                else if (barrier_hit is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag = 1;
                else if (barrier_hit is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag = 1;
                else if (barrier_hit is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag = 0;

                double OT = Math.Max(flag * option_type * (St - this.Strike), 0) + beta1 * cv;

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++)
                    {
                        double St = S;
                        double cv = 0;
                        List<double> path_St = new List<double>();

                        bool barrier_hit = false;
                        int flag = 0;

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double t = (i - 1)*dt;
                            double d1 = (Math.Log(St / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double delta = BlackScholesDelta(d1, this.IsCall);
                            double epsilon = RandNums[d,i];
                            double Stn = St * Math.Exp(nudt + sigsdt * epsilon);
                            cv = cv + delta * (Stn - St * erddt);
                            St = Stn;
                            path_St.Add(St);
                        }

                        if (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.DownAndIn)
                            barrier_hit = path_St.Where(x => x < this.BarrierLevel).Count() > 0;
                        else if (this.KnockType == KnockType.UpAndOut | this.KnockType == KnockType.UpAndIn)
                            barrier_hit = path_St.Where(x => x > this.BarrierLevel).Count() > 0;

                        if (barrier_hit is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag = 0;
                        else if (barrier_hit is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag = 1;
                        else if (barrier_hit is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag = 1;
                        else if (barrier_hit is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag = 0;

                        double OT = Math.Max(flag * option_type * (St - this.Strike), 0) + beta1 * cv;

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r*T);
            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M)*Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            double sum_OT = 0;
            double sum_OT2 = 0;

            int option_type = this.IsCall is true ? 1 : -1;

            for (int j = 0; j < M; j++)  // for each simulation
            {
                bool barrier_hit1 = false;
                bool barrier_hit2 = false;
                int flag1 = 0;
                int flag2 = 0;

                double St1 = S;
                double St2 = S;
                double cv1 = 0;
                double cv2 = 0;
                List<double> path_St1 = new List<double>();
                List<double> path_St2 = new List<double>();

                for (int i = 0; i < N; i++)  // for each time step
                {
                    double t = (i - 1) * dt;

                    double d1 = (Math.Log(St1 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    double d1_2 = (Math.Log(St2 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                    
                    double delta1 = BlackScholesDelta(d1, this.IsCall);
                    double delta2 = BlackScholesDelta(d1_2, this.IsCall);
                    
                    double epsilon = RandNums[j,i];

                    double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                    double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                    cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                    cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                    St1 = Stn1;
                    St2 = Stn2;
                    path_St1.Add(St1);
                    path_St2.Add(St2);
                }

                if (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.DownAndIn)
                {
                    barrier_hit1 = path_St1.Where(x => x < this.BarrierLevel).Count() > 0;
                    barrier_hit2 = path_St2.Where(x => x < this.BarrierLevel).Count() > 0;
                }
                else if (this.KnockType == KnockType.UpAndOut | this.KnockType == KnockType.UpAndIn)
                {
                    barrier_hit1 = path_St1.Where(x => x > this.BarrierLevel).Count() > 0;
                    barrier_hit2 = path_St2.Where(x => x > this.BarrierLevel).Count() > 0;
                }

                if (barrier_hit1 is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag1 = 0;
                else if (barrier_hit1 is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag1 = 1;
                else if (barrier_hit1 is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag1 = 1;
                else if (barrier_hit1 is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag1 = 0;

                if (barrier_hit2 is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag2 = 0;
                else if (barrier_hit2 is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                    flag2 = 1;
                else if (barrier_hit2 is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag2 = 1;
                else if (barrier_hit2 is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                    flag2 = 0;

                double OT = 0.5 * (Math.Max(flag1 * option_type * (St1 - this.Strike), 0) + beta1 * cv1 +
                                    Math.Max(flag2 * option_type * (St2 - this.Strike), 0) + beta1 + cv2);

                sum_OT = sum_OT + OT;
                sum_OT2 = sum_OT2 + OT*OT;
            }

            double option_value = sum_OT / M * Math.Exp(-r * T);
            double SD = Math.Sqrt( (sum_OT2 - sum_OT * sum_OT/M) * Math.Exp(-2*r*T)/(M-1));
            double SE = SD / Math.Sqrt(M);

            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }

        public List<double> GetPrice_Anti_ContVar_Multithread(double S, double T, double r, double div, double sig, int N, int M, double[,] RandNums)
        {
            double dt = T / N;

            double nudt = (r - div - 0.5 * Math.Pow(sig, 2))*dt;
            double sigsdt = sig*Math.Sqrt(dt);
            double erddt = Math.Exp((r - div) * dt);

            double beta1 = 0.5;

            List<double> sum_OT = new List<double>(new double[(int)M]);
            List<double> sum_OT2 = new List<double>(new double[(int)M]);

            int option_type = this.IsCall is true ? 1 : -1;

            // get number of processors in running computer
            int num_processors = System.Environment.ProcessorCount;

            // number of elements to generate for each processor
            int num_of_each = (int)(M / num_processors);

            // find remaining quantity; to be added to last Task
            int remainder = (int)M - num_of_each * num_processors;

            // set array of tasks by number of processors
            Task[] tasks = new Task[num_processors];

            for (int p = 0; p < num_processors; p++) {
                // to get by value instead of by reference
                int q = p;

                // starting index for each task
                int start_point = q * num_of_each;

                // ending index for each task
                int end_point = q * num_of_each + num_of_each;

                // set tasks array
                tasks[p] = new Task(() => {
                    // if the last task, add remainding quantity to num_of_each
                    if (q == num_processors - 1){
                        end_point = end_point + remainder;
                    }

                    // loop through num_of_each times
                    for (int d = start_point; d < end_point; d++) {
                        bool barrier_hit1 = false;
                        bool barrier_hit2 = false;
                        int flag1 = 0;
                        int flag2 = 0;

                        double St1 = S;
                        double St2 = S;
                        double cv1 = 0;
                        double cv2 = 0;
                        List<double> path_St1 = new List<double>();
                        List<double> path_St2 = new List<double>();

                        for (int i = 0; i < N; i++)  // for each time step
                        {
                            double t = (i - 1) * dt;

                            double d1 = (Math.Log(St1 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            double d1_2 = (Math.Log(St2 / this.Strike) + (r - div + 0.5 * Math.Pow(sig, 2)) * T) / (sig * Math.Sqrt(T));
                            
                            double delta1 = BlackScholesDelta(d1, this.IsCall);
                            double delta2 = BlackScholesDelta(d1_2, this.IsCall);
                            
                            double epsilon = RandNums[d,i];

                            double Stn1 = St1 * Math.Exp(nudt + sigsdt * epsilon);
                            double Stn2 = St2 * Math.Exp(nudt + sigsdt * -epsilon);

                            cv1 = cv1 + delta1 * (Stn1 - St1 * erddt);
                            cv2 = cv2 + delta2 * (Stn2 - St2 * erddt);
                            St1 = Stn1;
                            St2 = Stn2;
                            path_St1.Add(St1);
                            path_St2.Add(St2);
                        }

                        if (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.DownAndIn)
                        {
                            barrier_hit1 = path_St1.Where(x => x < this.BarrierLevel).Count() > 0;
                            barrier_hit2 = path_St2.Where(x => x < this.BarrierLevel).Count() > 0;
                        }
                        else if (this.KnockType == KnockType.UpAndOut | this.KnockType == KnockType.UpAndIn)
                        {
                            barrier_hit1 = path_St1.Where(x => x > this.BarrierLevel).Count() > 0;
                            barrier_hit2 = path_St2.Where(x => x > this.BarrierLevel).Count() > 0;
                        }

                        if (barrier_hit1 is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag1 = 0;
                        else if (barrier_hit1 is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag1 = 1;
                        else if (barrier_hit1 is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag1 = 1;
                        else if (barrier_hit1 is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag1 = 0;

                        if (barrier_hit2 is true & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag2 = 0;
                        else if (barrier_hit2 is false & (this.KnockType == KnockType.DownAndOut | this.KnockType == KnockType.UpAndOut))
                            flag2 = 1;
                        else if (barrier_hit2 is true & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag2 = 1;
                        else if (barrier_hit2 is false & (this.KnockType == KnockType.DownAndIn | this.KnockType == KnockType.UpAndIn))
                            flag2 = 0;

                        double OT = 0.5 * (Math.Max(flag1 * option_type * (St1 - this.Strike), 0) + beta1 * cv1 +
                                        Math.Max(flag2 * option_type * (St2 - this.Strike), 0) + beta1 + cv2);

                        sum_OT[d] = OT;
                        sum_OT2[d] = OT*OT;
                    }
                
                });
            }

            // call each task to run
            foreach (Task task in tasks) {
                task.Start();
                task.Wait();
            }

            double option_value = sum_OT.Sum() / M * Math.Exp(-r * T);

            double SD = Math.Sqrt( (sum_OT2.Sum() - sum_OT.Sum()*sum_OT.Sum()/M) * Math.Exp(-2 * r * T) / (M - 1) );
            double SE = SD / Math.Sqrt(M);

            // return option value and standard error in a list
            List<double> MyOutput = new List<double>();
            MyOutput.Add(option_value);
            MyOutput.Add(SE);
            return MyOutput;
        }
    }
}