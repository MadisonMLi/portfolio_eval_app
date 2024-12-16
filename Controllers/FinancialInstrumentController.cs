using Microsoft.AspNetCore.Mvc;
using System.Text.Json;

namespace Homework6;

[ApiController]
[Route("[controller]")]
public class FinancialInstrumentontroller : ControllerBase
{

    [HttpPost("/FinancialInstrument")]
    public ActionResult<FinancialInstrument> PostFinancialInstrument([FromBody] FinancialInstrument financialInstrument)
    {
        Console.WriteLine("Posted");
        FinanceContext db = new FinanceContext();
        Market mymarket = db.Markets.Single(x => x.Id == financialInstrument.marketid);

        modifydb.modifyfinancialinstrument(financialInstrument.Id, mymarket);
        return Ok(financialInstrument);
    }

    [HttpGet("/FinancialInstrument")]
    public ActionResult<FinancialInstrument> GetFinancialInstrument()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.FinancialInstruments.ToArray());
    }
}
