using Microsoft.AspNetCore.Mvc;
using System.Text.Json;

namespace Homework6;

[ApiController]
[Route("[controller]")]
public class Derivativesontroller : ControllerBase
{

    [HttpPost("/EuropeanPostData")]
    public ActionResult<European> PostEuropeanData([FromBody] European european)
    {
        Console.WriteLine("Posted");
        Console.WriteLine("before mod");

        FinanceContext db = new FinanceContext();

        Market mymarket = db.Markets.Single(x => x.Id == european.marketid);
        Underlying myunderlying = db.Underlyings.Single(x => x.Id == european.underlyingid);

        // set to correct datetime format
        european.expiration_date = DateTimeOffset.Parse(european.expiration_date.ToString()).UtcDateTime;
        Console.WriteLine(european.Is_Call);
        modifydb.modifyeuropean(mymarket, myunderlying, european.expiration_date, european.Strike, european.Is_Call);
        return Ok(european);
    }

    [HttpGet("/EuropeanGetData")]
    public ActionResult<European> GetEuropeanData()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.EuropeanOptions.ToArray());
    }

    [HttpPost("/DigitalPostData")]
    public ActionResult<European> PostDigitalData([FromBody] Digital digital)
    {
        Console.WriteLine("Posted");
        Console.WriteLine("before mod");

        FinanceContext db = new FinanceContext();

        Market mymarket = db.Markets.Single(x => x.Id == digital.marketid);
        Underlying myunderlying = db.Underlyings.Single(x => x.Id == digital.underlyingid);

        // set to correct datetime format
        digital.expiration_date = DateTimeOffset.Parse(digital.expiration_date.ToString()).UtcDateTime;

        modifydb.modifydigital(mymarket, myunderlying, digital.expiration_date, digital.Strike, digital.Is_Call, digital.Payout);
        return Ok(digital);
    }

    [HttpGet("/DigitalGetData")]
    public ActionResult<European> GetDigitalData()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.DigitalOptions.ToArray());
    }

    [HttpPost("/AsianPostData")]
    public ActionResult<European> PostAsianData([FromBody] Asian asian)
    {
        Console.WriteLine("Posted");
        Console.WriteLine("before mod");

        FinanceContext db = new FinanceContext();

        Market mymarket = db.Markets.Single(x => x.Id == asian.marketid);
        Underlying myunderlying = db.Underlyings.Single(x => x.Id == asian.underlyingid);

        // set to correct datetime format
        asian.expiration_date = DateTimeOffset.Parse(asian.expiration_date.ToString()).UtcDateTime;

        modifydb.modifyasian(mymarket, myunderlying, asian.expiration_date, asian.Strike, asian.Is_Call);
        return Ok(asian);
    }

    [HttpGet("/AsianGetData")]
    public ActionResult<European> GetAsianData()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.AsianOptions.ToArray());
    }

    [HttpPost("/RangePostData")]
    public ActionResult<European> PostRangeData([FromBody] Range range)
    {
        Console.WriteLine("Posted");
        Console.WriteLine("before mod");

        FinanceContext db = new FinanceContext();

        Market mymarket = db.Markets.Single(x => x.Id == range.marketid);
        Underlying myunderlying = db.Underlyings.Single(x => x.Id == range.underlyingid);

        // set to correct datetime format
        range.expiration_date = DateTimeOffset.Parse(range.expiration_date.ToString()).UtcDateTime;

        modifydb.modifyrange(mymarket, myunderlying, range.expiration_date);
        return Ok(range);
    }

    [HttpGet("/RangeGetData")]
    public ActionResult<European> GetRangeData()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.RangeOptions.ToArray());
    }

    [HttpPost("/LookbackPostData")]
    public ActionResult<European> PostLookbackData([FromBody] Lookback lookback)
    {
        Console.WriteLine("Posted");
        Console.WriteLine("before mod");

        FinanceContext db = new FinanceContext();

        Market mymarket = db.Markets.Single(x => x.Id == lookback.marketid);
        Underlying myunderlying = db.Underlyings.Single(x => x.Id == lookback.underlyingid);

        // set to correct datetime format
        lookback.expiration_date = DateTimeOffset.Parse(lookback.expiration_date.ToString()).UtcDateTime;

        modifydb.modifylookback(mymarket, myunderlying, lookback.expiration_date, lookback.Strike, lookback.Is_Call);
        return Ok(lookback);
    }

    [HttpGet("/LookbackGetData")]
    public ActionResult<European> GetLookbackData()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.LookbackOptions.ToArray());
    }

    [HttpPost("/BarrierPostData")]
    public ActionResult<European> PostBarrierData([FromBody] Barrier barrier)
    {
        Console.WriteLine("Posted");
        Console.WriteLine("before mod");

        FinanceContext db = new FinanceContext();

        Market mymarket = db.Markets.Single(x => x.Id == barrier.marketid);
        Underlying myunderlying = db.Underlyings.Single(x => x.Id == barrier.underlyingid);

        // set to correct datetime format
        barrier.expiration_date = DateTimeOffset.Parse(barrier.expiration_date.ToString()).UtcDateTime;

        modifydb.modifybarrier(mymarket, myunderlying, barrier.expiration_date, barrier.Strike, barrier.Is_Call, barrier.Barrier_Level, barrier.Knock_Type);
        return Ok(barrier);
    }

    [HttpGet("/BarrierGetData")]
    public ActionResult<European> GetBarrierData()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.BarrierOptions.ToArray());
    }
}
