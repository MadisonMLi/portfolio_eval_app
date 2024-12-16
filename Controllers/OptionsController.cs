using Microsoft.AspNetCore.Mvc;
using System.Text.Json;

namespace Homework6;

[ApiController]
[Route("[controller]")]
public class OptionsController : ControllerBase
{

    [HttpPost("/European")]
    public ActionResult<Options.UserOutputs> PostEuropean([FromBody] Options.UserInputs inputs)
    {
        Options.European makeoption = new Options.European(inputs.Strike, inputs.is_Call);
        double[,] rands = Options.GenRand(inputs.N, inputs.M);

        int antithetic_convert = Convert.ToInt32(inputs.antithetic);
        int multithread_convert = Convert.ToInt32(inputs.multithread);
        int cont_var_convert = Convert.ToInt32(inputs.cont_var);

        var func = makeoption.GetPrice;

        if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 0)
            func = makeoption.GetPrice_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar_Multithread;

        Console.WriteLine(inputs.Strike);
        Console.WriteLine(inputs.is_Call);
        Console.WriteLine("option price", makeoption.GetPrice(inputs.S, inputs.T, inputs.r, inputs.div, inputs.sig, inputs.N, inputs.M, rands));
        Console.WriteLine("{0}, {1}, {2}. {3}. {4}, {5}, {6}", inputs.S, inputs.T, inputs.r, inputs.div, inputs.sig, inputs.N, inputs.M);

        Dictionary<string, double> myoutput = Options.get_outputs(func,
                                                                  inputs.S,
                                                                  inputs.Strike,
                                                                  inputs.T,
                                                                  inputs.r,
                                                                  inputs.div,
                                                                  inputs.sig,
                                                                  inputs.N,
                                                                  inputs.M,
                                                                  rands,
                                                                  inputs.is_Call
                                                                //   inputs.OptionType,
                                                                //   inputs.Payout,
                                                                //   inputs.BarrierLevel,
                                                                //   inputs.knock,
                                                                //   inputs.VarReduction,
                                                                //   inputs.multithread
                                                                  );
        

        Options.UserOutputs user_outputs = new Options.UserOutputs(myoutput["Price"],
        myoutput["Standard Error"],
        myoutput["Delta"],
        myoutput["Gamma"],
        myoutput["Vega / 100"],
        myoutput["Theta / 365"],
        myoutput["Rho / 100"]
        );

        Console.WriteLine(JsonSerializer.Serialize(user_outputs));

        return Ok(user_outputs);
    }

    [HttpPost("/Asian")]
    public ActionResult<Options.UserOutputs> PostAsian([FromBody] Options.UserInputs inputs)
    {
        Options.Asian makeoption = new Options.Asian(inputs.Strike, inputs.is_Call);
        double[,] rands = Options.GenRand(inputs.N, inputs.M);

        int antithetic_convert = Convert.ToInt32(inputs.antithetic);
        int multithread_convert = Convert.ToInt32(inputs.multithread);
        int cont_var_convert = Convert.ToInt32(inputs.cont_var);

        var func = makeoption.GetPrice;

        if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 0)
            func = makeoption.GetPrice_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar_Multithread;

        Dictionary<string, double> myoutput = Options.get_outputs(func,
                                                                  inputs.S,
                                                                  inputs.Strike,
                                                                  inputs.T,
                                                                  inputs.r,
                                                                  inputs.div,
                                                                  inputs.sig,
                                                                  inputs.N,
                                                                  inputs.M,
                                                                  rands,
                                                                  inputs.is_Call);
                                                                //   inputs.OptionType,
                                                                //   inputs.Payout,
                                                                //   inputs.BarrierLevel,
                                                                //   inputs.knock,
                                                                //   inputs.VarReduction,
                                                                //   inputs.multithread

        Options.UserOutputs user_outputs = new Options.UserOutputs(myoutput["Price"],
        myoutput["Standard Error"],
        myoutput["Delta"],
        myoutput["Gamma"],
        myoutput["Vega / 100"],
        myoutput["Theta / 365"],
        myoutput["Rho / 100"]
        );

        Console.WriteLine(JsonSerializer.Serialize(user_outputs));

        return Ok(user_outputs);;
    }

    [HttpPost("/Digital")]
    public ActionResult<Options.UserOutputs> PostDigital([FromBody] Options.UserInputsDigital inputs)
    {
        Options.Digital makeoption = new Options.Digital(inputs.Strike, inputs.is_Call, inputs.Payout);
        double[,] rands = Options.GenRand(inputs.N, inputs.M);

        Console.WriteLine("Digital Option");
        Console.WriteLine("Payout {0}", inputs.Payout);

        int antithetic_convert = Convert.ToInt32(inputs.antithetic);
        int multithread_convert = Convert.ToInt32(inputs.multithread);
        int cont_var_convert = Convert.ToInt32(inputs.cont_var);

        var func = makeoption.GetPrice;

        if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 0)
            func = makeoption.GetPrice_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar_Multithread;

        Dictionary<string, double> myoutput = Options.get_outputs(func,
                                                                  inputs.S,
                                                                  inputs.Strike,
                                                                  inputs.T,
                                                                  inputs.r,
                                                                  inputs.div,
                                                                  inputs.sig,
                                                                  inputs.N,
                                                                  inputs.M,
                                                                  rands,
                                                                  inputs.is_Call);
                                                                //   inputs.OptionType,
                                                                //   inputs.Payout,
                                                                //   inputs.BarrierLevel,
                                                                //   inputs.knock,
                                                                //   inputs.VarReduction,
                                                                //   inputs.multithread

        Options.UserOutputs user_outputs = new Options.UserOutputs(myoutput["Price"],
        myoutput["Standard Error"],
        myoutput["Delta"],
        myoutput["Gamma"],
        myoutput["Vega / 100"],
        myoutput["Theta / 365"],
        myoutput["Rho / 100"]
        );

        Console.WriteLine(JsonSerializer.Serialize(user_outputs));

        return Ok(user_outputs);
    }

    [HttpPost("/Lookback")]
    public ActionResult<Options.UserOutputs> PostLookback([FromBody] Options.UserInputs inputs)
    {
        Options.Lookback makeoption = new Options.Lookback(inputs.Strike, inputs.is_Call);
        double[,] rands = Options.GenRand(inputs.N, inputs.M);

        int antithetic_convert = Convert.ToInt32(inputs.antithetic);
        int multithread_convert = Convert.ToInt32(inputs.multithread);
        int cont_var_convert = Convert.ToInt32(inputs.cont_var);

        var func = makeoption.GetPrice;

        if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 0)
            func = makeoption.GetPrice_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar_Multithread;

        Dictionary<string, double> myoutput = Options.get_outputs(func,
                                                                  inputs.S,
                                                                  inputs.Strike,
                                                                  inputs.T,
                                                                  inputs.r,
                                                                  inputs.div,
                                                                  inputs.sig,
                                                                  inputs.N,
                                                                  inputs.M,
                                                                  rands,
                                                                  inputs.is_Call);
                                                                //   inputs.OptionType,
                                                                //   inputs.Payout,
                                                                //   inputs.BarrierLevel,
                                                                //   inputs.knock,
                                                                //   inputs.VarReduction,
                                                                //   inputs.multithread

        Options.UserOutputs user_outputs = new Options.UserOutputs(myoutput["Price"],
        myoutput["Standard Error"],
        myoutput["Delta"],
        myoutput["Gamma"],
        myoutput["Vega / 100"],
        myoutput["Theta / 365"],
        myoutput["Rho / 100"]
        );

        Console.WriteLine(JsonSerializer.Serialize(user_outputs));

        return Ok(user_outputs);
    }

    [HttpPost("/Range")]
    public ActionResult<Options.UserOutputs> PostRange([FromBody] Options.UserInputs inputs)
    {
        Options.Range makeoption = new Options.Range();
        double[,] rands = Options.GenRand(inputs.N, inputs.M);

        int antithetic_convert = Convert.ToInt32(inputs.antithetic);
        int multithread_convert = Convert.ToInt32(inputs.multithread);
        int cont_var_convert = Convert.ToInt32(inputs.cont_var);

        var func = makeoption.GetPrice;

        if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 0)
            func = makeoption.GetPrice_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar_Multithread;

        Dictionary<string, double> myoutput = Options.get_outputs(func,
                                                                  inputs.S,
                                                                  inputs.Strike,
                                                                  inputs.T,
                                                                  inputs.r,
                                                                  inputs.div,
                                                                  inputs.sig,
                                                                  inputs.N,
                                                                  inputs.M,
                                                                  rands,
                                                                  inputs.is_Call);
                                                                //   inputs.OptionType,
                                                                //   inputs.Payout,
                                                                //   inputs.BarrierLevel,
                                                                //   inputs.knock,
                                                                //   inputs.VarReduction,
                                                                //   inputs.multithread

        Options.UserOutputs user_outputs = new Options.UserOutputs(myoutput["Price"],
        myoutput["Standard Error"],
        myoutput["Delta"],
        myoutput["Gamma"],
        myoutput["Vega / 100"],
        myoutput["Theta / 365"],
        myoutput["Rho / 100"]
        );

        Console.WriteLine(JsonSerializer.Serialize(user_outputs));

        return Ok(user_outputs);
    }

    [HttpPost("/Barrier")]
    public ActionResult<Options.UserOutputs> PostBarrier([FromBody] Options.UserInputsBarrier inputs)
    {
        Options.KnockType useKnock = Options.KnockType.DownAndOut;

        Console.WriteLine("Barrier Option");
        Console.WriteLine("Knock type integer {0}", inputs.knock);
        Console.WriteLine("Barrier level {0}", inputs.BarrierLevel);

        // decide on KnockType depending on integer value
        if (inputs.knock == 1)
            useKnock = Options.KnockType.UpAndOut;
        else if (inputs.knock == 2)
            useKnock = Options.KnockType.DownAndIn;
        else if (inputs.knock == 3)
            useKnock = Options.KnockType.UpAndIn;

        Options.Barrier makeoption = new Options.Barrier(inputs.Strike, inputs.is_Call, inputs.BarrierLevel, useKnock);
        double[,] rands = Options.GenRand(inputs.N, inputs.M);

        int antithetic_convert = Convert.ToInt32(inputs.antithetic);
        int multithread_convert = Convert.ToInt32(inputs.multithread);
        int cont_var_convert = Convert.ToInt32(inputs.cont_var);

        var func = makeoption.GetPrice;

        if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 0)
            func = makeoption.GetPrice_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 0)
            func = makeoption.GetPrice_Anti_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 0 & cont_var_convert == 1)
            func = makeoption.GetPrice_ContVar_Multithread;
        else if (multithread_convert == 0 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar;
        else if (multithread_convert == 1 & antithetic_convert == 1 & cont_var_convert == 1)
            func = makeoption.GetPrice_Anti_ContVar_Multithread;

        Dictionary<string, double> myoutput = Options.get_outputs(func,
                                                                  inputs.S,
                                                                  inputs.Strike,
                                                                  inputs.T,
                                                                  inputs.r,
                                                                  inputs.div,
                                                                  inputs.sig,
                                                                  inputs.N,
                                                                  inputs.M,
                                                                  rands,
                                                                  inputs.is_Call);
                                                                //   inputs.OptionType,
                                                                //   inputs.Payout,
                                                                //   inputs.BarrierLevel,
                                                                //   inputs.knock,
                                                                //   inputs.VarReduction,
                                                                //   inputs.multithread

        Options.UserOutputs user_outputs = new Options.UserOutputs(myoutput["Price"],
        myoutput["Standard Error"],
        myoutput["Delta"],
        myoutput["Gamma"],
        myoutput["Vega / 100"],
        myoutput["Theta / 365"],
        myoutput["Rho / 100"]
        );

        Console.WriteLine(JsonSerializer.Serialize(user_outputs));

        return Ok(user_outputs);
    }
}
