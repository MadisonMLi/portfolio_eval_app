using Microsoft.AspNetCore.Mvc;
using System.Text.Json;

namespace Homework6;

[ApiController]
[Route("[controller]")]
public class UnderlyingController : ControllerBase
{

    [HttpPost("/Underlying")]
    public ActionResult<Underlying> PostUnderlying([FromBody] Underlying underlying)
    {
        Console.WriteLine("Posted");

        FinanceContext db = new FinanceContext();

        Market mymarket = db.Markets.Single(x => x.Id == underlying.marketid);

        modifydb.modifyunderlyings(underlying.Id, underlying.marketid, underlying.Name, underlying.Symbol, mymarket);
        return Ok(underlying);

        // Console.WriteLine("Posted");
        // FinanceContext db = new FinanceContext();
        // Market mymarket = db.Markets.Single(x => x.Id == financialInstrument.marketid);

        // modifydb.modifyfinancialinstrument(financialInstrument.Id, mymarket);
        // return Ok(financialInstrument);
    }

    [HttpGet("/Underlying")]
    public ActionResult<Underlying> GetUnderlying()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.Underlyings.ToArray());
    }
}
