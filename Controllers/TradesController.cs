using Microsoft.AspNetCore.Mvc;
using System.Text.Json;

namespace Homework6;

[ApiController]
[Route("[controller]")]
public class TradeController : ControllerBase
{
    [HttpPost("/Trade")]
    public ActionResult<Trade> PostTrades([FromBody] Trade trade)
    {
        Console.WriteLine("Posted");

        FinanceContext db = new FinanceContext();

        modifydb.modifytrades(trade.Id, trade.quantity, trade.financialinstrumentid, trade.trade_price);
        return Ok(trade);
    }

    [HttpGet("/Trade")]
    public ActionResult<Trade> GetTrade()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.Trades.ToArray());
    }
}
