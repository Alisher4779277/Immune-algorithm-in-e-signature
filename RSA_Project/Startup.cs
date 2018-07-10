using Microsoft.Owin;
using Owin;

[assembly: OwinStartupAttribute(typeof(RSA_Project.Startup))]
namespace RSA_Project
{
    public partial class Startup
    {
        public void Configuration(IAppBuilder app)
        {
            ConfigureAuth(app);
        }
    }
}
